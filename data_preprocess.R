###This script is to implement preprocessing steps.

generateTimefrequence <- function(array.data, resolution, band.width,
                                  freq.list) {
  sample.size <- dim(array.data)
  mid.index <- floor(seq(band.width / 2 + 1, sample.size[1] - band.width / 2 - 1, length.out = resolution))
  index.matrix <- rbind(mid.index - band.width / 2, mid.index + band.width / 2)
  tf.array <- array(0, c(sample.size[2:3], length(freq.list), resolution))
  extreme.matrix <- matrix(0, sample.size[3], sample.size[2])
  for (tri in 1:sample.size[3]) {
    for (chan in 1:sample.size[2]) {
        tf.array[chan, tri, , ] <- generateSingleperiodogram(index.matrix, array.data[, chan, tri] / sd(array.data[, chan, tri]), name.bands, freq.list)
        extreme.matrix[tri, chan] <- min(tf.array[chan, tri, , ])
    }
  }
  return(tf.array)
}

generateSingleperiodogram <- function(index.matrix, single.signals, name.bands,
                                      freq.list) {
  time.freq <- c()
  for (j in 1:ncol(index.matrix)) {
    single.truncated <- single.signals[index.matrix[1, j]:index.matrix[2, j]]
    periodogram <- (abs(fft(single.truncated)/sqrt(length(single.truncated))) ^ 2)
    power.vec <- c()
    for (band in name.bands) {
      power.vec <- c(power.vec, mean(periodogram[floor(unlist(freq.list[band]) * length(single.truncated) / 1000)]))
    }
    time.freq <- cbind(time.freq, power.vec)
  }
  return(time.freq)
}
