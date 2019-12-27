# echo "setup PmAlg PmAlg-00-00-01 in /workfs/bes/leizh/workarea-6.6.5.p01/Analysis/Physics"

if ( $?CMTROOT == 0 ) then
  setenv CMTROOT /afs/ihep.ac.cn/bes3/offline/ExternalLib/SLC6/contrib/CMT/v1r25
endif
source ${CMTROOT}/mgr/setup.csh
set cmtPmAlgtempfile=`${CMTROOT}/mgr/cmt -quiet build temporary_name`
if $status != 0 then
  set cmtPmAlgtempfile=/tmp/cmt.$$
endif
${CMTROOT}/mgr/cmt setup -csh -pack=PmAlg -version=PmAlg-00-00-01 -path=/workfs/bes/leizh/workarea-6.6.5.p01/Analysis/Physics  -no_cleanup $* >${cmtPmAlgtempfile}
if ( $status != 0 ) then
  echo "${CMTROOT}/mgr/cmt setup -csh -pack=PmAlg -version=PmAlg-00-00-01 -path=/workfs/bes/leizh/workarea-6.6.5.p01/Analysis/Physics  -no_cleanup $* >${cmtPmAlgtempfile}"
  set cmtsetupstatus=2
  /bin/rm -f ${cmtPmAlgtempfile}
  unset cmtPmAlgtempfile
  exit $cmtsetupstatus
endif
set cmtsetupstatus=0
source ${cmtPmAlgtempfile}
if ( $status != 0 ) then
  set cmtsetupstatus=2
endif
/bin/rm -f ${cmtPmAlgtempfile}
unset cmtPmAlgtempfile
exit $cmtsetupstatus

