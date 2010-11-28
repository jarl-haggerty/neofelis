(ns com.curious.neofelis.extend-orfs
  (:use clojure.java.io
        clojure.java.shell
        com.curious.neofelis.utils)
  (:require [clojure.string :as string]
            [clojure.xml :as xml]
            [clojure.zip :as zip]
            [clojure.contrib.zip-filter.xml :as zf]))

(defn- parse-blast-xml
  "This function will parse the xml that blast outputs as a result of blasting the possible extensions of genes found by genemark.  The result of this function will be
   a map that uses the original locations of the genes found by genemark as keys, and values which are vectors that contain the information to put into a excel file
   should an extension prove to be better than the original"
  [blast-xml]
  (apply merge
         (for [iteration (zf/xml-> blast-xml :BlastOutput_iterations :Iteration)
               :when (seq (zf/xml-> iteration :Iteration_hits :Hit :Hit_hsps :Hsp))]
           (let [hsps (seq (zf/xml-> iteration :Iteration_hits :Hit :Hit_hsps :Hsp))
                 [original-location new-location] (rest (re-find #"_(\d+-\d+)_(\d+-\d+)"
                                                                 (first (zf/xml-> iteration :Iteration_query-def zf/text))))
                 [best-hit best-hsp] (let [result (apply min-key #(Double/valueOf (first (zf/xml-> % :Hsp_evalue zf/text))) hsps)]
                                       [(-> result zip/up zip/up) result])
                 [title organism] (let [raw-hit-def (first (zf/xml-> best-hit :Hit_def zf/text))
                                        gt-index (.indexOf raw-hit-def ">")
                                        hit-def (.substring raw-hit-def 0 (if (> gt-index 0)
                                                                            gt-index
                                                                            (count raw-hit-def)))]
                                    (if (> (.indexOf hit-def "[") 0)
                                      [(.trim (.substring hit-def 0 (.indexOf hit-def "[")))
                                       (.trim (.replace (.substring hit-def (inc (.indexOf hit-def "[")))
                                                        "]" ""))]
                                      [(.trim hit-def) (.trim hit-def)]))]
             {original-location [new-location
                                 (count (zf/xml-> iteration :Iteration_hits :Hit))
                                 (first (zf/xml-> best-hsp :Hsp_bit-score zf/text))
                                 (first (zf/xml-> best-hsp :Hsp_evalue zf/text))
                                 (first (zf/xml-> best-hsp :Hsp_identity zf/text))
                                 (first (zf/xml-> best-hsp :Hsp_align-len zf/text))
                                 (first (zf/xml-> best-hit :Hit_id zf/text))
                                 title
                                 organism]}))))

(defn extend-genes
  "Accepts the location of a blast installation(blast), a database(database), a file name(query), and a name for the genome in the query(name).
   This function will extract the start, stop, and e-value for every gene from the original genemark prediction from the spreadsheet that was generated.
   Then for every gene this function will search from just before the start of the gene, to either the beginning of the sequence or the end of a gene
   just before it, for a start codon.  If at any time a stop codon is encountered the search stops.  The protein sequence of each extension found this
   was is then written to a canidates.fas file which is then blasted."
  [blast database query name e-value]
  (when-not (.exists (file "extended-blasts"))
    (.mkdir (file "extended-blasts")))
  (when-not (.exists (file "extended-blast-spreadsheets"))
      (.mkdir (file "extended-blast-spreadsheets")))
  (let [;this will read in the relevent data from the spreadsheet made
        ;by the initial gene prediction and blast.  forward-orfs will
        ;be a series of vectors, [start stop e-value], and reverse
        ;orfs will be series of vectors [stop start e-value].
        spreadsheet (file (str "initial-blast-spreadsheets/" name ".blastp.xls"))
        [forward-orfs reverse-orfs] (with-open [input (reader spreadsheet)]
                                      (let [data (->> input line-seq rest (map #(string/split % #"\t")))
                                            locations-evalues (map #(conj (vec (string/split (nth % 1) #"-")) (nth % 4)) data)
                                            orfs (map #(vector (Integer/valueOf (first %))
                                                               (Integer/valueOf (second %))
                                                               (try (Double/valueOf (last %))
                                                                    (catch NumberFormatException _ Double/POSITIVE_INFINITY)))
                                                      locations-evalues)]
                                        [(doall (map #(vector (-> % first dec) (second %) (last %))
                                                     (filter #(< (first %) (second %)) orfs)))
                                         (doall (map #(vector (-> % second dec) (first %) (last %))
                                                     (filter #(> (first %) (second %)) orfs)))]))
        ;forward-stops and reverse-stops represent the boundaries of
        ;how far back a search for the new start of a gene can go.
        forward-stops (map #(- (second %) 2) forward-orfs)
        reverse-stops (map #(+ (first %) 1) reverse-orfs)
        genome (with-open [input (reader query)]
                 (->> input
                      line-seq
                      (map string/upper-case)
                      (filter #(re-find #"^[ATCG]+$" %))
                      string/join))
        ;forward-extensions is a list of all the alternative starts of
        ;a gene.  The first item in the sequence is the original gene,
        ;repesented as [start stop e-value.  The rest of the sequence
        ;is alternative starts found by looking at the codons
        ;before the start of the original sequence and after the
        ;latest forward-stop before the start of this sequence, or the
        ;beginning of the genome if there is no such forward stop.
        forward-extensions (for [orf forward-orfs]
                             (apply vector orf
                                    (for [start (range (- (first orf) 3)
                                                       (let [lower-stops (filter #(< % (first orf)) forward-stops)]
                                                         (apply max 0 lower-stops))
                                                       -3)
                                          :while (not (contains? stop-codons (.substring genome start (+ start 3))))
                                          :when (contains? start-codons (.substring genome start (+ start 3)))]
                                      start)))
        ;reverse-extensions is found with the same method as that for
        ;the forward-extensions, except the algorithm is performed on
        ;the reverse complement genome with the reverse-stops.
        reverse-extensions (for [orf reverse-orfs]
                             (apply vector orf
                                    (for [start (range (second orf)
                                                       (let [upper-stops (filter #(> % (second orf)) reverse-stops)]
                                                         (apply min (count genome) upper-stops))
                                                       3)
                                          :while (not (contains? stop-codons (reverse-complement (.substring genome start (+ start 3)))))
                                          :when (contains? start-codons (reverse-complement (.substring genome start (+ start 3))))]
                                      start)))]
    ;This form writes all the possible extensions found into
    ;canidates.fas and the next form blasts this file if there is no
    ;xml file in the extended-blasts directory for this genome.
    (with-open [output (writer "canidates.fas")]
      (doseq [extension forward-extensions]
        (let [original-zero-location (first extension)
              original-fasta-start (inc (get original-zero-location 0))
              original (assoc original-zero-location 0 original-fasta-start)]
          (doseq [start (rest extension)]
            (.write output (str ">canidate_" (string/join "-" (take 2 original)) "_" start "-" (second original) "\n"))
            (.write output (string/join "\n" (map #(apply str %) (partition-all 50 (-> genome (.substring start (second original)) translate)))))
            (.newLine output))))
      (doseq [extension reverse-extensions]
        (let [original-zero-location (first extension)
              original-fasta-start (inc (get original-zero-location 0))
              original (assoc original-zero-location 0 original-fasta-start)]
          (doseq [start (rest extension)]
            (.write output (str ">canidate_" (string/join "-" (->> original (take 2) reverse)) "_" start "-" (first original) "\n"))
            (.write output (string/join "\n" (map #(apply str %) (partition-all 50 (-> genome (.substring (first original) start) reverse-complement translate)))))
            (.newLine output)))))
    (when-not (.exists (file (str "extended-blasts/" name ".blastp.xml")))
      (println "Blasting")
      (with-open [output (writer (str "extended-blasts/" name ".blastp.xml"))]
        (let [result (:out (sh (str blast "/bin/blastp")
                               "-db" (str blast "/db/" database)
                               "-num_threads" (str (.availableProcessors (Runtime/getRuntime)))
                               "-evalue" (str e-value)
                               "-outfmt" "5"
                               "-query" "canidates.fas"))]
          (.write output (.substring result (.indexOf result "<BlastOutput>"))))))
    (println "Blasted")
    (let [blast-xml (-> (str "extended-blasts/" name ".blastp.xml") xml/parse zip/xml-zip)
          good-extensions (parse-blast-xml blast-xml)
          original-excel (with-open [input (reader (str "initial-blast-spreadsheets/" name ".blastp.xls"))]
                           (doall (->> input
                                       line-seq
                                       rest
                                       (map #(string/split % #"\t")))))]
      (with-open [output (writer (str "extended-blast-spreadsheets/" name ".blastp.xls"))]
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
        (doseq [line original-excel]
          (if (and (contains? good-extensions (second line))
                   (< (try (Double/valueOf (nth (get good-extensions (second line)) 3))
                           (catch NumberFormatException _ Double/POSITIVE_INFINITY))
                      (try (Double/valueOf (nth line 4))
                           (catch NumberFormatException _ Double/POSITIVE_INFINITY))))
            (->> (second line)
                 (get good-extensions)
                 (cons (first line))
                 (string/join "\t")
                 (.write output))
            (.write output (string/join "\t" line)))
          (.write output "\n"))))))
