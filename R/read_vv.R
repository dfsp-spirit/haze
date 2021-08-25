#' @title Read vv binary file.
#'
#' @description Read matrix-like data from vv files. This is a custom format from the cpp_geodesic repo, designed to allow fast reading of vector-of-vectors data. The format does NOT require that all inner vectors have the same length, so it is NOT limited to matrices. The format is designed for storing graphs as adjacency lists.
#'
#' @param filepath string. Full path to the input vv file.
#'
#' @param datatype one of `integer()` or `float()`, suitable for the given file.
#'
#' @return list of vectors, the data.
#'
#' @export
read.vv <- function(filepath, datatype=integer()) {
  MAGIC_FILE_TYPE_NUMBER_1 = 42;
  MAGIC_FILE_TYPE_NUMBER_2 = 13;
  endian = "big";

  fh = file(filepath, "rb");
  on.exit({ close(fh) }, add=TRUE);

  magic1 = readBin(fh, integer(), n = 1, size = 4, endian = endian);
  magic2 = readBin(fh, integer(), n = 1, size = 4, endian = endian);
  if (magic1 != MAGIC_FILE_TYPE_NUMBER_1) {
    stop(sprintf("Magic number mismatch (%d != %d).", magic1, MAGIC_FILE_TYPE_NUMBER_1)); # nocov
  }
  if (magic2 != MAGIC_FILE_TYPE_NUMBER_2) {
    stop(sprintf("Magic number mismatch (%d != %d).", magic2, MAGIC_FILE_TYPE_NUMBER_2)); # nocov
  }

  data = list();
  num_vectors = readBin(fh, integer(), n = 1, size = 4, endian = endian);
  for(vec_idx in seq.int(num_vectors)) {
    this_vec_len = readBin(fh, integer(), n = 1, size = 4, endian = endian);
    data[[length(data)+1L]] = readBin(fh, datatype, n = this_vec_len, size = 4, endian = endian);
  }
  return(data);
}

