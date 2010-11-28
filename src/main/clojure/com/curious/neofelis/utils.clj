(ns com.curious.neofelis.utils
  (:require [clojure.contrib.zip-filter.xml :as zf]
            [clojure.zip :as zip]
            [clojure.string :as string]))

(def
 ^{:doc "Start codonss"}
 start-codons #{"ATG" "GTG" "TTG"})
(def
 ^{:doc "Stop codonss"}
 stop-codons #{"TGA" "TAA" "TAG"})

(def
 ^{:doc "A map with contains codons as keys and the cooresponding amino acids as values"}
 translation-map
     (apply hash-map
            (concat (interleave ["TTT" "TTC"] (repeat 2 "F"))
                    (interleave ["TTA" "TTG" "CTT" "CTC" "CTA" "CTG"] (repeat 6 "L"))
                    (interleave ["TCT" "TCC" "TCA" "TCG"] (repeat 4 "S"))
                    (interleave ["TAT" "TAC"] (repeat 2 "Y"))
                    (interleave ["TGA" "TAA" "TAG"] (repeat 3 "*"))
                    (interleave ["TGT" "TGC"] (repeat 2 "C"))
                    (interleave ["TGG"] (repeat 1 "W"))
                    (interleave ["CCA" "CCC" "CCG" "CCT"] (repeat 4 "P"))
                    (interleave ["CAC" "CAT"] (repeat 2 "H"))
                    (interleave ["CAA" "CAG"] (repeat 2 "Q"))
                    (interleave ["CGA" "CGC" "CGG" "CGT" "AGA" "AGG"] (repeat 6 "R"))
                    (interleave ["ATT" "ATC" "ATA"] (repeat 3 "I"))
                    (interleave ["ATG"] (repeat 1 "M"))
                    (interleave ["ACA", "ACC", "ACG", "ACT"] (repeat 4 "T"))
                    (interleave ["AAT", "AAC"] (repeat 2 "N"))
                    (interleave ["AAA", "AAG"] (repeat 2 "K"))
                    (interleave ["AGT", "AGC"] (repeat 2 "S"))
                    (interleave ["GTA", "GTT", "GTG", "GTC"] (repeat 4 "V"))
                    (interleave ["GCA", "GCT", "GCG", "GCC"] (repeat 4 "A"))
                    (interleave ["GAT", "GAC"] (repeat 2 "D"))
                    (interleave ["GAA" "GAG"] (repeat 2 "E"))
                    (interleave ["GGA", "GGT", "GGG", "GGC"] (repeat 4 "G")))))

(defn translate
  "Translates a nucleotide sequence(a string) into an amino acid sequence(another string)."
  [input]
  (->> input
       (partition 3)
       (map #(apply str %))
       (map #(get translation-map %))
       (apply str)))

(defn reverse-complement [input]
  "Returns the reverse complement of a nucleotide sequence."
  (let [complement-map {\A \T \T \A \G \C \C \G}
        result (string/escape input #(get complement-map %))]
    (string/reverse result)))
