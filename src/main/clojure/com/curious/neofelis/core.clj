(ns com.curious.neofelis.core
  (:import java.io.File)
  (:require [clojure.string :as string])
  (:use clojure.java.io
        clojure.contrib.command-line
        com.curious.neofelis.phage-pipeline
        com.curious.neofelis.extend-orfs)
  (:gen-class))

(defn expand-directories 
  "This function takes a list of files and returns a list with all directories removed and the files in those directories added."
  [input]
  (filter (fn [x] (not (.isHidden x)))
          (flatten (map (fn [x] (if (.isDirectory x) (seq (.listFiles x)) x))
                        (map file input)))))

(def
 ^{:doc "Extensions recognized by the annotation miner."}
 extensions #{".seq" ".faa" ".fna" ".fas" ".fasta" ".genome"})

(defn -main [& args]
  "This is the main method of the annotation miner.

     The following sequence of events occurs for each genome:

     If the settings contain a genemark matrix then genemark is called on the
     genome using that matrix.  The results from this call are then used as input
     to a blastp search using the evalue defined in settings.  Finally, an excel
     spreadsheet is created to summarize the results of the blast search."
  (with-command-line args
    "annotation-miner"
    [[evalue "pass an evalue to use when running the blast search on the genemark predictions" "0.1"]
     [matrix "pass a matrix to use when predicting genes with genemark" nil]
     [heuristically? "perform the genemark predictions heuristically"]
     [cutoffs "cutoffs" "0"]
     [min-length "pass the minimum length of an intergenic gene" "100"]
     [max-e "max-e" "1000"]
     [genemark "location of genemark installation" "~/genemark"]
     [blast "location of blast installation" "~/blast"]
     [database "database to run blast with" "nr"]
     [prefix "prefix" ""]
     queries]

    (doseq [genome (expand-directories queries)]
      (let [name (let [extension-start (.lastIndexOf (.getName genome) ".")]
                   (if (and (< -1 extension-start) (contains? extensions (.substring (.getName genome) extension-start)))
                     (do (.substring (.getName genome) 0 extension-start))
                     (.getName genome)))]
        (when (or heuristically? matrix)
          (let [pipeline-arguments [matrix
                                    (.getAbsolutePath genome)
                                    name
                                    (str name prefix)
                                    database
                                    (.replace blast "~" (System/getProperty "user.home"))
                                    (.replace genemark "~" (System/getProperty "user.home"))
                                    (Double/valueOf evalue)]]
            (apply phage-pipeline (if heuristically? (rest pipeline-arguments) pipeline-arguments))))
        (extend-genes (.replace blast "~" (System/getProperty "user.home")) database (.getAbsolutePath genome) name evalue)
          (comment
          (intergenic (:min-length settings)
                      (:blastall settings)
                      (:database settings)
                      (:max-e settings)))))))


