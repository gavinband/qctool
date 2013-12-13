bgen.read_snp_id_data
<- function(header) {
    require(bitops)
    con = header$connection
    endian = header$endian
    data = list()
    data$number_of_samples = readBin(con, integer(), size = 4, 
        n = 1, endian = endian)
    if (length(data$number_of_samples) == 0) {
        return(NULL)
    }
    id_storage_size = readBin(con, integer(), size = 1, n = 1, 
        endian = endian, signed = FALSE)
    SNPID_size = readBin(con, integer(), size = 1, n = 1, 
        endian = endian, signed = FALSE)
    data$SNPID = readChar(con, nchars = SNPID_size)
    if (id_storage_size > SNPID_size) {
        readChar(con, nchar = id_storage_size - SNPID_size)
    }
    rsid_size = readBin(con, integer(), size = 1, n = 1, 
        endian = endian, signed = FALSE)
    data$rsid = readChar(con, nchars = rsid_size)
    if (id_storage_size > rsid_size) {
        readChar(con, nchar = id_storage_size - rsid_size)
    }
    chromosome = readBin(con, integer(), size = 1, n = 1, 
        endian = endian, signed = FALSE)
    if (chromosome < 23) {
        data$chromosome = sprintf("%02d", chromosome)
    }
    else if (chromosome == 23) {
        data$chromosome = "0X"
    }
    else if (chromosome == 24) {
        data$chromosome = "0Y"
    }
    else if (chromosome == 253) {
        data$chromosome = "XY"
    }
    else if (chromosome == 254) {
        data$chromosome = "MT"
    }
    else if (chromosome == 255) {
        data$chromosome = NA
    }
    else {
        stop(sprintf("Unknown chromosome %d", chromosome))
    }
    data$position = readBin(con, integer(), size = 4, n = 1, 
        endian = endian)
    if (bitAnd(header$flags, 2)) {
        alleleA_size = readBin(con, integer(), size = 1, 
            n = 1, endian = endian, signed = FALSE)
        data$alleleA = readChar(con, nchars = alleleA_size)
        alleleB_size = readBin(con, integer(), size = 1, 
            n = 1, endian = endian, signed = FALSE)
        data$alleleB = readChar(con, nchars = alleleB_size)
    }
    else {
        data$alleleA = readChar(con, nchars = 1)
        data$alleleB = readChar(con, nchars = 1)
    }
    return(data)
}
