#' @title Read vv binary file.
#'
#' @description Read matrix-like data from vv files. This is a custom format from the cpp_geodesic repo, designed to allow fast reading of vector-of-vectors data. The format does NOT require that all inner vectors have the same length, so it is NOT limited to matrices. The format is designed for storing graphs as adjacency lists.
#'
#' @param filepath string. Full path to the input vv file.
#'
#' @return list of vectors, the data. The vv files may can store double or int, which is encoded in the file header and used accordingly.
#'
#' @export
read.vv <- function(filepath) {
  MAGIC_FILE_TYPE_NUMBER = 42;
  DATA_TYPE_INT32 = 13;
  DATA_TYPE_FLOAT32 = 14;
  endian = "big";

  fh = file(filepath, "rb");
  on.exit({ close(fh) }, add=TRUE);

  magic = readBin(fh, integer(), n = 1, size = 4, endian = endian);
  data_type_code = readBin(fh, integer(), n = 1, size = 4, endian = endian);
  if (magic != MAGIC_FILE_TYPE_NUMBER) {
    stop(sprintf("Magic number mismatch (%d != %d).", magic, MAGIC_FILE_TYPE_NUMBER)); # nocov
  }
  if (!(data_type_code %in% c(DATA_TYPE_INT32, DATA_TYPE_FLOAT32))) {
    stop(sprintf("Invalid data type code '%d', supported are '13' for int32 and '14' for float32.\n", data_type_code)); # nocov
  }

  if(data_type_code == DATA_TYPE_INT32) {
    datatype = integer();
  } else {
    datatype = numeric();
  }

  data = list();
  num_vectors = readBin(fh, integer(), n = 1, size = 4, endian = endian);
  for(vec_idx in seq.int(num_vectors)) {
    this_vec_len = readBin(fh, integer(), n = 1, size = 4, endian = endian);
    data[[length(data)+1L]] = readBin(fh, datatype, n = this_vec_len, size = 4, endian = endian);
  }
  return(data);
}

