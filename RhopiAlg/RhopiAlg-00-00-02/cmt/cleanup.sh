# echo "cleanup omega RhopiAlg-00-00-02 in /afs/ihep.ac.cn/users/l/leizh"

if test "${CMTROOT}" = ""; then
  CMTROOT=/afs/ihep.ac.cn/bes3/offline/ExternalLib/SLC6/contrib/CMT/v1r25; export CMTROOT
fi
. ${CMTROOT}/mgr/setup.sh
cmtomegatempfile=`${CMTROOT}/mgr/cmt -quiet build temporary_name`
if test ! $? = 0 ; then cmtomegatempfile=/tmp/cmt.$$; fi
${CMTROOT}/mgr/cmt cleanup -sh -pack=omega -version=RhopiAlg-00-00-02 -path=/afs/ihep.ac.cn/users/l/leizh  $* >${cmtomegatempfile}
if test $? != 0 ; then
  echo >&2 "${CMTROOT}/mgr/cmt cleanup -sh -pack=omega -version=RhopiAlg-00-00-02 -path=/afs/ihep.ac.cn/users/l/leizh  $* >${cmtomegatempfile}"
  cmtcleanupstatus=2
  /bin/rm -f ${cmtomegatempfile}
  unset cmtomegatempfile
  return $cmtcleanupstatus
fi
cmtcleanupstatus=0
. ${cmtomegatempfile}
if test $? != 0 ; then
  cmtcleanupstatus=2
fi
/bin/rm -f ${cmtomegatempfile}
unset cmtomegatempfile
return $cmtcleanupstatus

