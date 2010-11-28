(ns com.curious.neofelis.intergenic
  (:use clojure.java.io
        com.curious.neofelis.utils)
  (:require [clojure.string :as string]))

(defn- parse-blast-xml
  [blast-xml]
  (for [iteration (zf/xml-> blast-xml :BlastOutput_iterations :Iteration) :when (seq (zf/xml-> iteration :Iteration_hits :Hit :Hit_hsps :Hsp))]
    (when-let [hsps (seq (zf/xml-> iteration :Iteration_hits :Hit :Hit_hsps :Hsp))]
      ;(println hsps)
      (let [;_ (println (zf/xml-> iteration :Iteration_query-def zf/text))
            [query-id query-location] (string/split (first (zf/xml-> iteration :Iteration_query-def zf/text)) #":")
            ;_ (println query-id query-location)
            [best-hit best-hsp] (let [result (apply min-key #(Double/valueOf (first (zf/xml-> % :Hsp_evalue zf/text))) hsps)]
                                  [(-> result zip/up zip/up) result])
            ;_ (println (-> best-hit zip/node :tag) (-> best-hsp zip/node :tag))
            [title organism] (let [raw-hit-def (first (zf/xml-> best-hit :Hit_def zf/text))
                                   gt-index (.indexOf raw-hit-def ">")
                                   hit-def (.substring raw-hit-def 0 (if (> gt-index 0) gt-index (count raw-hit-def)))]
                               ;(println hit-def)
                               (if (> (.indexOf hit-def "[") 0)
                                 [(.trim (.substring hit-def 0 (.indexOf hit-def "[")))
                                  (.trim (.replace (.substring hit-def (inc (.indexOf hit-def "[")))
                                                   "]" ""))]
                                 [(.trim hit-def) (.trim hit-def)]))
                                        ;_ (println title "|||" organism)
            ]
        [query-id
         query-location
         (count (zf/xml-> iteration :Iteration_hits :Hit))
         (first (zf/xml-> best-hsp :Hsp_bit-score zf/text))
         (first (zf/xml-> best-hsp :Hsp_evalue zf/text))
         (first (zf/xml-> best-hsp :Hsp_identity zf/text))
         (first (zf/xml-> best-hsp :Hsp_align-len zf/text))
         (first (zf/xml-> best-hit :Hit_id zf/text))
         title
         organism]))))

(defn remove-from-interval [remaining deletions min-length]
  (if-let [deletion (first deletions)]
    (recur (apply concat (for [section remaining]
                           (filter #(> (second %) (first %))
                                   [(if (< (first deletion) (second section))
                                      [(first section) (first deletion)]
                                      [(first section) (second section)])
                                    (if (> (second deletion) (first section))
                                      [(second deletion) (second section)]
                                      [(first section) (second section)])])))
           (rest deletions)
           min-length)
    remaining))

(defn find-intergenics [min-length blast database max-e name query]
  (let [genes (with-open [input (reader (str "extended-blast-spreadsheets/" name ".blastp.xls"))]
                (doall (for [line (rest (line-seq input))]
                         (string/split line #"\t"))))
        genome (with-open [input (reader query)]
                 (->> input
                      line-seq
                      (map string/upper-case)
                      (filter #(re-find #"^[ATCG]+$" %))
                      string/join))
        [forward-filled-regions reverse-filled-regions] (apply map #(->> %& (keep identity) vec)
                                                               (for [gene genes]
                                                                 (let [pieces (map #(Integer/valueOf %)
                                                                                   (string/split (second gene) #"-"))]
                                                                   (if (< (first pieces) (second pieces))
                                                                     [[(dec (first pieces)) (second pieces)] nil]
                                                                     [nil [(dec (second pieces)) (first pieces)]]))))
        [forward-empty-regions reverse-empty-regions] (map remove-from-interval
                                                           [[[0 (count genome)]] [[0 (count genome)]]]
                                                           [forward-filled-regions reverse-filled-regions]
                                                           [min-length min-length])
        forward-canidates (filter #(> (- (second %) (first %)) min-length)
                                  (for [region forward-empty-regions
                                        frame (range 3)
                                        start (range (- (first region) frame) (- (second region) 1))
                                        :when (contains? start-codons (.substring genome start (+ start 3)))]
                                    [start (first (for [stop (range (+ start 3) (- (second region) 1))
                                                        :when (contains? stop-codons (.substring genome stop (+ stop 3)))]
                                                    stop))]))
        reverse-canidates (filter #(> (- (second %) (first %)) min-length)
                                  (for [region reverse-empty-regions
                                        frame (range 3)
                                        start (range (+ (second region) frame) (- (first region) 2))
                                        :when (contains? start-codons (reverse-complement (.substring genome (- start 3) start)))]
                                    [start (first (for [stop (range (- start 3) (- (second region) 2))
                                                        :when (contains? stop-codons (.substring genome (- stop 3) stop))]
                                                    stop))]))
        ]
    (with-open [output (writer (str "intergenic-canidates/" name ".intergenics.fas"))]
      (doseq [canidate forward-canidates]
        (.write output (str ">canidate_" (inc (first canidate)) "-" (second canidate) "\n"))
        (.write output (partition-all 50 (.substring genome (first canidate) (second canidate)))))
      (doseq [canidate reverse-canidates]
        (.write output (str ">canidates_" (second canidate) "-" (inc (first canidate)) "\n"))
        (.write output (partition-all 50 (.substring genome (first canidate) (second canidate))))))
    (when-not (.exists (file (str "intergenic-blasts/" name ".blastp.xml")))
      (println "Blasting")
      (with-open [output (writer (str "intergenic-blasts/" name ".blastp.xml"))]
        (let [result (:out (sh (str blast "/bin/blastp")
                               "-db" (str blast "/db/" database)
                               "-num_threads" (str (.availableProcessors (Runtime/getRuntime)))
                               "-evalue" (str e-value)
                               "-outfmt" "5"
                               "-query" (str "intergenic-canidates/" name ".intergenics.fas")))]
          (.write output (.substring result (.indexOf result "<BlastOutput>")))))
      (println "Blasted"))
    (let [blast-xml (xml/parse (str "intergenic-blasts/" name ".blastp.xml"))
          canidate-hits (filter #(< (Double/valueOf  (nth % 4)) max-e) (parse-blast-xml blast-xml))]
      (with-open [output (writer (file (str "intergenic-blast-spreadsheets/" name ".blastp.xls")))]
           (.write output (str (string/join "\t" ["query-id"
                                            "query-location"
                                            "total-hits"
                                            "best-bit-score"
                                            "best-e-value"
                                            "best-identity"
                                            "Hsp-align-len"
                                            "best-hit-gi"
                                            "definition"
                                            "organism"]) "\n"))
           (doseq [gene canidate-hits]
             (.write output (str (string/join "\t" gene) "\n")))))))
