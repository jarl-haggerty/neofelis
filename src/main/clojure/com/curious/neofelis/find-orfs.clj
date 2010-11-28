(ns com.curious.neofelis.find-orfs)

(def file-filter
     (proxy [Object FileFilter] []
       (accept [file] (.isHidden file))))

(defn find-orfs [prefix & raw-queries]
  (let [queries (flatten (map (fn [x] (if (.isDirectory x) (seq (.listFiles x file-filter)) x))
                              (map as-file queries)))
        orfs-directory (as-file "orfs")]
    (if (not (.isDirectory orfs-directory))
      (.mkdir orfs-directory))))
