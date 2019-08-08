# echo "setup omega RhopiAlg-00-00-02 in /afs/ihep.ac.cn/users/l/leizh"

if ( $?CMTROOT == 0 ) then
  setenv CMTROOT /afs/ihep.ac.cn/bes3/offline/ExternalLib/SLC6/contrib/CMT/v1r25
endif
source ${CMTROOT}/mgr/setup.csh
set cmtomegatempfile=`${CMTROOT}/mgr/cmt -quiet build temporary_name`
if $status != 0 then
  set cmtomegatempfile=/tmp/cmt.$$
endif
${CMTROOT}/mgr/cmt setup -csh -pack=omega -version=RhopiAlg-00-00-02 -path=/afs/ihep.ac.cn/users/l/leizh  -no_cleanup $* >${cmtomegatempfile}
if ( $status != 0 ) then
  echo "${CMTROOT}/mgr/cmt setup -csh -pack=omega -version=RhopiAlg-00-00-02 -path=/afs/ihep.ac.cn/users/l/leizh  -no_cleanup $* >${cmtomegatempfile}"
  set cmtsetupstatus=2
  /bin/rm -f ${cmtomegatempfile}
  unset cmtomegatempfile
  exit $cmtsetupstatus
endif
set cmtsetupstatus=0
source ${cmtomegatempfile}
if ( $status != 0 ) then
  set cmtsetupstatus=2
endif
/bin/rm -f ${cmtomegatempfile}
unset cmtomegatempfile
exit $cmtsetupstatus

