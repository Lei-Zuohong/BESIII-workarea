# echo "setup omega RhopiAlg-00-00-02 in /afs/ihep.ac.cn/users/l/leizh"

if test "${CMTROOT}" = ""; then
  CMTROOT=/afs/ihep.ac.cn/bes3/offline/ExternalLib/SLC6/contrib/CMT/v1r25; export CMTROOT
fi
. ${CMTROOT}/mgr/setup.sh
cmtomegatempfile=`${CMTROOT}/mgr/cmt -quiet build temporary_name`
if test ! $? = 0 ; then cmtomegatempfile=/tmp/cmt.$$; fi
${CMTROOT}/mgr/cmt setup -sh -pack=omega -version=RhopiAlg-00-00-02 -path=/afs/ihep.ac.cn/users/l/leizh  -no_cleanup $* >${cmtomegatempfile}
if test $? != 0 ; then
  echo >&2 "${CMTROOT}/mgr/cmt setup -sh -pack=omega -version=RhopiAlg-00-00-02 -path=/afs/ihep.ac.cn/users/l/leizh  -no_cleanup $* >${cmtomegatempfile}"
  cmtsetupstatus=2
  /bin/rm -f ${cmtomegatempfile}
  unset cmtomegatempfile
  return $cmtsetupstatus
fi
cmtsetupstatus=0
. ${cmtomegatempfile}
if test $? != 0 ; then
  cmtsetupstatus=2
fi
/bin/rm -f ${cmtomegatempfile}
unset cmtomegatempfile
return $cmtsetupstatus

