## code to prepare `nirsdata`
f = system.file("data-raw", package = "nirsls")|>
  list.files(pattern = ".txt")

dfs = lapply(f,
             read_portamon)

dfs[[1]] <- set_starttime(dfs[[1]], 14,2)

plot(dfs[[1]])

nirsdata <- lapply(f, read_rawnirs) |>
  data.table::rbindlist(
    idcol = NULL,
    use.names = TRUE
  )

usethis::use_data(nirsdata, overwrite = TRUE)

## code to prepare `pdata`
pdata = fs::path_package(
  "data-raw",
  "participantdata.xlsx",
  package = "nirsls") |>
  readxl::read_xlsx()

usethis::use_data(pdata, overwrite = TRUE)
