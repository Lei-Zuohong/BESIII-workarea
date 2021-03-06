# echo "setup OmegaAlg OmegaAlg-00-00-01 in /besfs/users/leizh/besiii/workarea665/Analysis/Physics"

if test "${CMTROOT}" = ""; then
  CMTROOT=/cvmfs/bes3.ihep.ac.cn/bes3sw/ExternalLib/SLC6/contrib/CMT/v1r25; export CMTROOT
fi
. ${CMTROOT}/mgr/setup.sh
cmtOmegaAlgtempfile=`${CMTROOT}/mgr/cmt -quiet build temporary_name`
if test ! $? = 0 ; then cmtOmegaAlgtempfile=/tmp/cmt.$$; fi
${CMTROOT}/mgr/cmt setup -sh -pack=OmegaAlg -version=OmegaAlg-00-00-01 -path=/besfs/users/leizh/besiii/workarea665/Analysis/Physics  -no_cleanup $* >${cmtOmegaAlgtempfile}
if test $? != 0 ; then
  echo >&2 "${CMTROOT}/mgr/cmt setup -sh -pack=OmegaAlg -version=OmegaAlg-00-00-01 -path=/besfs/users/leizh/besiii/workarea665/Analysis/Physics  -no_cleanup $* >${cmtOmegaAlgtempfile}"
  cmtsetupstatus=2
  /bin/rm -f ${cmtOmegaAlgtempfile}
  unset cmtOmegaAlgtempfile
  return $cmtsetupstatus
fi
cmtsetupstatus=0
. ${cmtOmegaAlgtempfile}
if test $? != 0 ; then
  cmtsetupstatus=2
fi
/bin/rm -f ${cmtOmegaAlgtempfile}
unset cmtOmegaAlgtempfile
return $cmtsetupstatus

