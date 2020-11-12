# echo "setup OmegaAlg OmegaAlg-00-00-01 in /besfs/users/leizh/besiii/workarea665/Analysis/Physics"

if ( $?CMTROOT == 0 ) then
  setenv CMTROOT /cvmfs/bes3.ihep.ac.cn/bes3sw/ExternalLib/SLC6/contrib/CMT/v1r25
endif
source ${CMTROOT}/mgr/setup.csh
set cmtOmegaAlgtempfile=`${CMTROOT}/mgr/cmt -quiet build temporary_name`
if $status != 0 then
  set cmtOmegaAlgtempfile=/tmp/cmt.$$
endif
${CMTROOT}/mgr/cmt setup -csh -pack=OmegaAlg -version=OmegaAlg-00-00-01 -path=/besfs/users/leizh/besiii/workarea665/Analysis/Physics  -no_cleanup $* >${cmtOmegaAlgtempfile}
if ( $status != 0 ) then
  echo "${CMTROOT}/mgr/cmt setup -csh -pack=OmegaAlg -version=OmegaAlg-00-00-01 -path=/besfs/users/leizh/besiii/workarea665/Analysis/Physics  -no_cleanup $* >${cmtOmegaAlgtempfile}"
  set cmtsetupstatus=2
  /bin/rm -f ${cmtOmegaAlgtempfile}
  unset cmtOmegaAlgtempfile
  exit $cmtsetupstatus
endif
set cmtsetupstatus=0
source ${cmtOmegaAlgtempfile}
if ( $status != 0 ) then
  set cmtsetupstatus=2
endif
/bin/rm -f ${cmtOmegaAlgtempfile}
unset cmtOmegaAlgtempfile
exit $cmtsetupstatus

