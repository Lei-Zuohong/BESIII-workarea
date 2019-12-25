#-- start of make_header -----------------

#====================================
#  Document PmAlg_check_install_runtime
#
#   Generated Tue Dec 10 17:45:00 2019  by leizh
#
#====================================

include ${CMTROOT}/src/Makefile.core

ifdef tag
CMTEXTRATAGS = $(tag)
else
tag       = $(CMTCONFIG)
endif

cmt_PmAlg_check_install_runtime_has_no_target_tag = 1

#--------------------------------------------------------

ifdef cmt_PmAlg_check_install_runtime_has_target_tag

tags      = $(tag),$(CMTEXTRATAGS),target_PmAlg_check_install_runtime

PmAlg_tag = $(tag)

#cmt_local_tagfile_PmAlg_check_install_runtime = $(PmAlg_tag)_PmAlg_check_install_runtime.make
cmt_local_tagfile_PmAlg_check_install_runtime = $(bin)$(PmAlg_tag)_PmAlg_check_install_runtime.make

else

tags      = $(tag),$(CMTEXTRATAGS)

PmAlg_tag = $(tag)

#cmt_local_tagfile_PmAlg_check_install_runtime = $(PmAlg_tag).make
cmt_local_tagfile_PmAlg_check_install_runtime = $(bin)$(PmAlg_tag).make

endif

include $(cmt_local_tagfile_PmAlg_check_install_runtime)
#-include $(cmt_local_tagfile_PmAlg_check_install_runtime)

ifdef cmt_PmAlg_check_install_runtime_has_target_tag

cmt_final_setup_PmAlg_check_install_runtime = $(bin)setup_PmAlg_check_install_runtime.make
cmt_dependencies_in_PmAlg_check_install_runtime = $(bin)dependencies_PmAlg_check_install_runtime.in
#cmt_final_setup_PmAlg_check_install_runtime = $(bin)PmAlg_PmAlg_check_install_runtimesetup.make
cmt_local_PmAlg_check_install_runtime_makefile = $(bin)PmAlg_check_install_runtime.make

else

cmt_final_setup_PmAlg_check_install_runtime = $(bin)setup.make
cmt_dependencies_in_PmAlg_check_install_runtime = $(bin)dependencies.in
#cmt_final_setup_PmAlg_check_install_runtime = $(bin)PmAlgsetup.make
cmt_local_PmAlg_check_install_runtime_makefile = $(bin)PmAlg_check_install_runtime.make

endif

#cmt_final_setup = $(bin)setup.make
#cmt_final_setup = $(bin)PmAlgsetup.make

#PmAlg_check_install_runtime :: ;

dirs ::
	@if test ! -r requirements ; then echo "No requirements file" ; fi; \
	  if test ! -d $(bin) ; then $(mkdir) -p $(bin) ; fi

javadirs ::
	@if test ! -d $(javabin) ; then $(mkdir) -p $(javabin) ; fi

srcdirs ::
	@if test ! -d $(src) ; then $(mkdir) -p $(src) ; fi

help ::
	$(echo) 'PmAlg_check_install_runtime'

binobj = 
ifdef STRUCTURED_OUTPUT
binobj = PmAlg_check_install_runtime/
#PmAlg_check_install_runtime::
#	@if test ! -d $(bin)$(binobj) ; then $(mkdir) -p $(bin)$(binobj) ; fi
#	$(echo) "STRUCTURED_OUTPUT="$(bin)$(binobj)
endif

${CMTROOT}/src/Makefile.core : ;
ifdef use_requirements
$(use_requirements) : ;
endif

#-- end of make_header ------------------
#-- start of cmt_action_runner_header ---------------

ifdef ONCE
PmAlg_check_install_runtime_once = 1
endif

ifdef PmAlg_check_install_runtime_once

PmAlg_check_install_runtimeactionstamp = $(bin)PmAlg_check_install_runtime.actionstamp
#PmAlg_check_install_runtimeactionstamp = PmAlg_check_install_runtime.actionstamp

PmAlg_check_install_runtime :: $(PmAlg_check_install_runtimeactionstamp)
	$(echo) "PmAlg_check_install_runtime ok"
#	@echo PmAlg_check_install_runtime ok

#$(PmAlg_check_install_runtimeactionstamp) :: $(PmAlg_check_install_runtime_dependencies)
$(PmAlg_check_install_runtimeactionstamp) ::
	$(silent) /afs/ihep.ac.cn/bes3/offline/Boss/6.6.5.p01/BesPolicy/BesPolicy-01-05-05/cmt/bes_check_installations.sh -files= -s=../share *.txt   -installdir=/workfs/bes/leizh/workarea-6.6.5.p01/InstallArea/share
	$(silent) cat /dev/null > $(PmAlg_check_install_runtimeactionstamp)
#	@echo ok > $(PmAlg_check_install_runtimeactionstamp)

PmAlg_check_install_runtimeclean ::
	$(cleanup_silent) /bin/rm -f $(PmAlg_check_install_runtimeactionstamp)

else

#PmAlg_check_install_runtime :: $(PmAlg_check_install_runtime_dependencies)
PmAlg_check_install_runtime ::
	$(silent) /afs/ihep.ac.cn/bes3/offline/Boss/6.6.5.p01/BesPolicy/BesPolicy-01-05-05/cmt/bes_check_installations.sh -files= -s=../share *.txt   -installdir=/workfs/bes/leizh/workarea-6.6.5.p01/InstallArea/share

endif

install ::
uninstall ::

#-- end of cmt_action_runner_header -----------------
#-- start of cleanup_header --------------

clean :: PmAlg_check_install_runtimeclean ;
#	@cd .

ifndef PEDANTIC
.DEFAULT::
	$(echo) "(PmAlg_check_install_runtime.make) $@: No rule for such target" >&2
else
.DEFAULT::
	$(error PEDANTIC: $@: No rule for such target)
endif

PmAlg_check_install_runtimeclean ::
#-- end of cleanup_header ---------------