# echo "setup PppmpzAlg PppmpzAlg-00-00-01 in /workfs/bes/leizh/workarea-6.6.5.p01/Analysis/Physics"

if test "${CMTROOT}" = ""; then
  CMTROOT=/cvmfs/bes3.ihep.ac.cn/bes3sw/ExternalLib/SLC6/contrib/CMT/v1r25; export CMTROOT
fi
. ${CMTROOT}/mgr/setup.sh
cmtPppmpzAlgtempfile=`${CMTROOT}/mgr/cmt -quiet build temporary_name`
if test ! $? = 0 ; then cmtPppmpzAlgtempfile=/tmp/cmt.$$; fi
${CMTROOT}/mgr/cmt setup -sh -pack=PppmpzAlg -version=PppmpzAlg-00-00-01 -path=/workfs/bes/leizh/workarea-6.6.5.p01/Analysis/Physics  -no_cleanup $* >${cmtPppmpzAlgtempfile}
if test $? != 0 ; then
  echo >&2 "${CMTROOT}/mgr/cmt setup -sh -pack=PppmpzAlg -version=PppmpzAlg-00-00-01 -path=/workfs/bes/leizh/workarea-6.6.5.p01/Analysis/Physics  -no_cleanup $* >${cmtPppmpzAlgtempfile}"
  cmtsetupstatus=2
  /bin/rm -f ${cmtPppmpzAlgtempfile}
  unset cmtPppmpzAlgtempfile
  return $cmtsetupstatus
fi
cmtsetupstatus=0
. ${cmtPppmpzAlgtempfile}
if test $? != 0 ; then
  cmtsetupstatus=2
fi
/bin/rm -f ${cmtPppmpzAlgtempfile}
unset cmtPppmpzAlgtempfile
return $cmtsetupstatus

