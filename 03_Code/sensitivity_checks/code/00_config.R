options(stringsAsFactors=FALSE, warn=1)
args0 <- grep('^--file=', commandArgs(), value=TRUE)
script0 <- if(length(args0)) sub('^--file=', '', args0[1]) else file.path(getwd(),'00_config.R')
ROOT <- Sys.getenv('VS_REPAIR_ROOT', unset=dirname(dirname(normalizePath(script0))))
BASE <- Sys.getenv('VS_BASELINE_ROOT', unset=normalizePath(file.path(ROOT,'..','..'),mustWork=FALSE))
DATA <- Sys.getenv('VS_SERVER_ROOT', unset=file.path(BASE,'input'))
RDATA <- file.path(DATA,'RData')
SC_INPUT <- Sys.getenv('VS_SC_QS', unset=file.path(BASE,'input','sc_nc_raw.qs'))
TPM_INPUT <- Sys.getenv('VS_TPM_XLSX', unset=file.path(BASE,'input','TPM.xlsx'))
EXTERNAL_ROOT <- Sys.getenv('VS_EXTERNAL_ROOT', unset=file.path(BASE,'input','external'))
OUT <- file.path(ROOT,'results'); LOG <- file.path(ROOT,'logs')
dir.create(OUT,recursive=TRUE,showWarnings=FALSE);dir.create(LOG,recursive=TRUE,showWarnings=FALSE)
tasklib <- file.path(ROOT,'R_library');dir.create(tasklib,recursive=TRUE,showWarnings=FALSE)
compatlib <- Sys.getenv('VS_SC_COMPAT_LIBRARY',unset='')
.libPaths(unique(c(tasklib,compatlib[nzchar(compatlib)],.libPaths())))
write_result <- function(x,name) write.csv(x,file.path(OUT,paste0(name,'.csv')),row.names=FALSE,na='')
load_env <- function(name){e<-new.env();load(file.path(RDATA,name),e);e}
start_log <- function(name){sink(file.path(LOG,paste0(name,'.log')),split=TRUE);cat('Started',format(Sys.time(),'%Y-%m-%d %H:%M:%S %Z'),'\n');cat(R.version.string,'\n')}
end_log <- function(){cat('Completed',format(Sys.time(),'%Y-%m-%d %H:%M:%S %Z'),'\n');print(sessionInfo());sink()}
