# echo "setup PppmpzAlg PppmpzAlg-00-00-01 in /besfs/users/leizh/besiii/workarea705/Analysis/Physics"

if ( $?CMTROOT == 0 ) then
  setenv CMTROOT /cvmfs/bes3.ihep.ac.cn/bes3sw/ExternalLib/SLC6/contrib/CMT/v1r25
endif
source ${CMTROOT}/mgr/setup.csh
set cmtPppmpzAlgtempfile=`${CMTROOT}/mgr/cmt -quiet build temporary_name`
if $status != 0 then
  set cmtPppmpzAlgtempfile=/tmp/cmt.$$
endif
${CMTROOT}/mgr/cmt setup -csh -pack=PppmpzAlg -version=PppmpzAlg-00-00-01 -path=/besfs/users/leizh/besiii/workarea705/Analysis/Physics  -no_cleanup $* >${cmtPppmpzAlgtempfile}
if ( $status != 0 ) then
  echo "${CMTROOT}/mgr/cmt setup -csh -pack=PppmpzAlg -version=PppmpzAlg-00-00-01 -path=/besfs/users/leizh/besiii/workarea705/Analysis/Physics  -no_cleanup $* >${cmtPppmpzAlgtempfile}"
  set cmtsetupstatus=2
  /bin/rm -f ${cmtPppmpzAlgtempfile}
  unset cmtPppmpzAlgtempfile
  exit $cmtsetupstatus
endif
set cmtsetupstatus=0
source ${cmtPppmpzAlgtempfile}
if ( $status != 0 ) then
  set cmtsetupstatus=2
endif
/bin/rm -f ${cmtPppmpzAlgtempfile}
unset cmtPppmpzAlgtempfile
exit $cmtsetupstatus

