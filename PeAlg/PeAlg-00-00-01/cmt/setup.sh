# echo "setup PmAlg PmAlg-00-00-01 in /workfs/bes/leizh/workarea-6.6.5.p01/Analysis/Physics"

if test "${CMTROOT}" = ""; then
  CMTROOT=/afs/ihep.ac.cn/bes3/offline/ExternalLib/SLC6/contrib/CMT/v1r25; export CMTROOT
fi
. ${CMTROOT}/mgr/setup.sh
cmtPmAlgtempfile=`${CMTROOT}/mgr/cmt -quiet build temporary_name`
if test ! $? = 0 ; then cmtPmAlgtempfile=/tmp/cmt.$$; fi
${CMTROOT}/mgr/cmt setup -sh -pack=PmAlg -version=PmAlg-00-00-01 -path=/workfs/bes/leizh/workarea-6.6.5.p01/Analysis/Physics  -no_cleanup $* >${cmtPmAlgtempfile}
if test $? != 0 ; then
  echo >&2 "${CMTROOT}/mgr/cmt setup -sh -pack=PmAlg -version=PmAlg-00-00-01 -path=/workfs/bes/leizh/workarea-6.6.5.p01/Analysis/Physics  -no_cleanup $* >${cmtPmAlgtempfile}"
  cmtsetupstatus=2
  /bin/rm -f ${cmtPmAlgtempfile}
  unset cmtPmAlgtempfile
  return $cmtsetupstatus
fi
cmtsetupstatus=0
. ${cmtPmAlgtempfile}
if test $? != 0 ; then
  cmtsetupstatus=2
fi
/bin/rm -f ${cmtPmAlgtempfile}
unset cmtPmAlgtempfile
return $cmtsetupstatus

