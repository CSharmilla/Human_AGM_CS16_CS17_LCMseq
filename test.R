
random_text = function(n, n_words = 1) {
  sapply(1:n, function(i) {
    w = replicate(sample(n_words, 1), 
                  paste0(sample(letters, sample(4:10, 1)), collapse = ""))
    if(n_words > 1) {
      paste(w, collapse = " ")
    } else {
      w
    }
  })
}


split = sample(letters[1:5], 100, replace = TRUE)
sentences = lapply(unique(split), function(x) {
  random_text(3, 8)
})
names(sentences) = unique(split)

mat = matrix(rnorm(100*10), nrow = 100)

Heatmap(mat, name = "mat", cluster_rows = FALSE,
        right_annotation = rowAnnotation(
          textbox = anno_textbox(
            list("a" = 1:10, "b" = 20:50), 
            sentences[c("a", "b")])
        )
)
