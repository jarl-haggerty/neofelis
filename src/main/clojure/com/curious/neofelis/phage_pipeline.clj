(ns com.curious.neofelis.phage-pipeline
  (:use clojure.java.io
        clojure.set
        clojure.java.shell)
  (:require [clojure.string :as string]
            [clojure.xml :as xml]
            [clojure.zip :as zip]
            [clojure.contrib.zip-filter.xml :as zf]))

(defn modify-fasta-header [input-file output-file label]
  (with-open [input (reader input-file)
              output (writer output-file)]
    (doseq [line (line-seq input)]
      (if-let [matches (next (re-find #"(>orf_\d+).*, (\d+ - \d+)" line))]
        (.write output 
                (-> (str (first matches) ":" (second matches) "\n")
                  (string/replace "orf" label)
                  (string/replace " " "")))
        (.write output (str line "\n"))))))

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

(defn phage-pipeline
  ([matrix query name prefix database blast genemark evalue]
     "This function will use the specified genemark executeable and matrix to predict genes within the specified query.  
      Then the blastall executable will be used to perform a blastp search on the predicting genes using the specified evalue.
      Finally, an excel speadsheet will be created to summarize the results of the process."
     (let [blast-directory (file "initial-blasts")
           excel-directory (file "initial-blast-spreadsheets")]
       (when-not (.isDirectory blast-directory)
         (.mkdir blast-directory))
       (when-not (.isDirectory excel-directory)
         (.mkdir excel-directory))
       
       (sh (str genemark "/gm") "-opq" "-m" matrix query)
       (modify-fasta-header (str query ".orf") (str prefix ".faa") name)
       (sh "rm" (str query ".orf"))
       
       (when-not (.exists (file (str "initial-blasts/" prefix ".blastp.xml")))
         (with-open [output (writer (str "initial-blasts/" prefix ".blastp.xml"))]
           (let [result (:out (sh (str blast "/bin/blastp")
                                  "-db" (str blast "/db/" database)
                                  "-num_threads" (str (.availableProcessors (Runtime/getRuntime)))
                                  "-evalue" (str evalue)
                                  "-outfmt" "5"
                                  "-query" (str prefix ".faa")))]
             (.write output (.substring result (.indexOf result "<BlastOutput>"))))))

       (let [blast-xml (zip/xml-zip (xml/parse (str "initial-blasts/" prefix ".blastp.xml")))
             genes (set (with-open [fasta (reader (str prefix ".faa"))]
                          (doall (for [line (line-seq fasta) :when (re-find #">" line)]
                                   (.substring line 1)))))
             good-genes (set (zf/xml-> blast-xml
                                       :BlastOutput_iterations
                                       :Iteration
                                       #(> (count (zf/xml-> % :Iteration_hits :Hit)) 0)
                                       :Iteration_query-def
                                       zf/text))]
         (with-open [output (writer (file (str "initial-blast-spreadsheets/" prefix ".blastp.xls")))]
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
           (doseq [gene (parse-blast-xml blast-xml)]
             (.write output (str (string/join "\t" gene) "\n")))
           (doseq [gene (difference genes good-genes)]
             (.write output (str (string/join "\t" [(.substring gene 0 (.indexOf gene ":"))
                                             (.substring gene (inc (.indexOf gene ":")))
                                             "0"
                                             "None"
                                             "Infinite"
                                             ">0.1"
                                             "None"
                                             "None"
                                             "None"
                                             "No hits <0.1"
                                             "None"]) "\n")))))))
([query name prefix database blast genemark evalue]
   (let [gc-content (last (re-find #"GC% = (\d+)"
                                   (:out (sh (str genemark "/gc") query))))]
     (phage-pipeline (str genemark "/heuristic_mat/heu_11_" gc-content ".mat") query name prefix database blast genemark evalue))))
