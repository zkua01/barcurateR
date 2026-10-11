# Install and manage barcurateR reference databases

Functions to locate, download, and manage cached reference databases.

- `rb_db_dir()`: Returns the path to the local cache directory.

- `rb_db_path()`: Returns the full path to a specific database file in
  the cache.

- `rb_db_available()`: Checks if a database file exists locally.

- `rb_db_url()`: Returns the configured download URL for the full
  database.

- `rb_install_db()`: Downloads the full database from Zenodo (or a
  custom URL) into the local cache, with retry logic and checksum
  verification.

- `rb_yzfishdb_release()`: Returns metadata for the official YZFishDB
  Zenodo release.

- `rb_default_yzfishdb_url()`: Returns the default Zenodo API URL for
  YZFishDB.

## Usage

``` r
rb_yzfishdb_release()

rb_default_yzfishdb_url()

rb_db_url()

rb_db_dir(create = TRUE)

rb_db_path(filename = "YZFishDB.db")

rb_db_available(path = rb_db_path())

rb_install_db(
  url = rb_db_url(),
  destfile = rb_db_path(),
  overwrite = FALSE,
  md5 = NULL,
  sha256 = NULL,
  min_bytes = 1,
  quiet = FALSE,
  timeout = 3600,
  retries = 3
)
```

## Arguments

- create:

  Logical. If `TRUE`, create the cache directory if it doesn't exist.

- filename:

  The name of the database file (default `"YZFishDB.db"`).

- path:

  Path to check for database availability.

- url:

  URL to download the database from.

- destfile:

  Destination file path for the download.

- overwrite:

  Logical. Overwrite existing file?

- md5, sha256:

  Expected checksums for verification.

- min_bytes:

  Minimum expected file size in bytes.

- quiet:

  Logical. Suppress download progress messages?

- timeout:

  Download timeout in seconds (default 3600).

- retries:

  Number of download retry attempts.

## Value

- `rb_db_dir()` and `rb_db_path()` return a character string path.

- `rb_db_available()` returns a logical value.

- `rb_install_db()` returns the normalized path to the downloaded file.

- `rb_yzfishdb_release()` returns a data frame of release metadata.

## Details

The cache directory is determined by:

1.  The `BARCURATER_CACHE_DIR` environment variable.

2.  The default R user data directory
    (`tools::R_user_dir("barcurateR", "data")`).

The download URL is determined by:

1.  The `BARCURATER_DB_URL` environment variable.

2.  The `barcurateR.db_url` global option.

3.  The default YZFishDB Zenodo URL.

## See also

Other database management:
[`rb_connect()`](https://zkua01.github.io/barcurateR/reference/rb_connect.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Check if the database is already installed
rb_db_available()

# Get the path where the database will be installed
rb_db_path()

# Install the YZFishDB database (downloads from Zenodo)
rb_install_db()
} # }
```
