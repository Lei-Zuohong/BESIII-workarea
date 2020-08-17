#-- start of make_header -----------------

#====================================
#  Document PppmpzAlg_check_install_runtime
#
#   Generated Tue Aug 11 07:11:31 2020  by leizh
#
#====================================

include ${CMTROOT}/src/Makefile.core

ifdef tag
CMTEXTRATAGS = $(tag)
else
tag       = $(CMTCONFIG)
endif

cmt_PppmpzAlg_check_install_runtime_has_no_target_tag = 1

#--------------------------------------------------------

ifdef cmt_PppmpzAlg_check_install_runtime_has_target_tag

tags      = $(tag),$(CMTEXTRATAGS),target_PppmpzAlg_check_install_runtime

PppmpzAlg_tag = $(tag)

#cmt_local_tagfile_PppmpzAlg_check_install_runtime = $(PppmpzAlg_tag)_PppmpzAlg_check_install_runtime.make
cmt_local_tagfile_PppmpzAlg_check_install_runtime = $(bin)$(PppmpzAlg_tag)_PppmpzAlg_check_install_runtime.make

else

tags      = $(tag),$(CMTEXTRATAGS)

PppmpzAlg_tag = $(tag)

#cmt_local_tagfile_PppmpzAlg_check_install_runtime = $(PppmpzAlg_tag).make
cmt_local_tagfile_PppmpzAlg_check_install_runtime = $(bin)$(PppmpzAlg_tag).make

endif

include $(cmt_local_tagfile_PppmpzAlg_check_install_runtime)
#-include $(cmt_local_tagfile_PppmpzAlg_check_install_runtime)

ifdef cmt_PppmpzAlg_check_install_runtime_has_target_tag

cmt_final_setup_PppmpzAlg_check_install_runtime = $(bin)setup_PppmpzAlg_check_install_runtime.make
cmt_dependencies_in_PppmpzAlg_check_install_runtime = $(bin)dependencies_PppmpzAlg_check_install_runtime.in
#cmt_final_setup_PppmpzAlg_check_install_runtime = $(bin)PppmpzAlg_PppmpzAlg_check_install_runtimesetup.make
cmt_local_PppmpzAlg_check_install_runtime_makefile = $(bin)PppmpzAlg_check_install_runtime.make

else

cmt_final_setup_PppmpzAlg_check_install_runtime = $(bin)setup.make
cmt_dependencies_in_PppmpzAlg_check_install_runtime = $(bin)dependencies.in
#cmt_final_setup_PppmpzAlg_check_install_runtime = $(bin)PppmpzAlgsetup.make
cmt_local_PppmpzAlg_check_install_runtime_makefile = $(bin)PppmpzAlg_check_install_runtime.make

endif

#cmt_final_setup = $(bin)setup.make
#cmt_final_setup = $(bin)PppmpzAlgsetup.make

#PppmpzAlg_check_install_runtime :: ;

dirs ::
	@if test ! -r requirements ; then echo "No requirements file" ; fi; \
	  if test ! -d $(bin) ; then $(mkdir) -p $(bin) ; fi

javadirs ::
	@if test ! -d $(javabin) ; then $(mkdir) -p $(javabin) ; fi

srcdirs ::
	@if test ! -d $(src) ; then $(mkdir) -p $(src) ; fi

help ::
	$(echo) 'PppmpzAlg_check_install_runtime'

binobj = 
ifdef STRUCTURED_OUTPUT
binobj = PppmpzAlg_check_install_runtime/
#PppmpzAlg_check_install_runtime::
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
PppmpzAlg_check_install_runtime_once = 1
endif

ifdef PppmpzAlg_check_install_runtime_once

PppmpzAlg_check_install_runtimeactionstamp = $(bin)PppmpzAlg_check_install_runtime.actionstamp
#PppmpzAlg_check_install_runtimeactionstamp = PppmpzAlg_check_install_runtime.actionstamp

PppmpzAlg_check_install_runtime :: $(PppmpzAlg_check_install_runtimeactionstamp)
	$(echo) "PppmpzAlg_check_install_runtime ok"
#	@echo PppmpzAlg_check_install_runtime ok

#$(PppmpzAlg_check_install_runtimeactionstamp) :: $(PppmpzAlg_check_install_runtime_dependencies)
$(PppmpzAlg_check_install_runtimeactionstamp) ::
	$(silent) /cvmfs/bes3.ihep.ac.cn/bes3sw/Boss/6.6.5.p01/BesPolicy/BesPolicy-01-05-05/cmt/bes_check_installations.sh -files= -s=../share *.txt   -installdir=/workfs/bes/leizh/workarea-6.6.5.p01/InstallArea/share
	$(silent) cat /dev/null > $(PppmpzAlg_check_install_runtimeactionstamp)
#	@echo ok > $(PppmpzAlg_check_install_runtimeactionstamp)

PppmpzAlg_check_install_runtimeclean ::
	$(cleanup_silent) /bin/rm -f $(PppmpzAlg_check_install_runtimeactionstamp)

else

#PppmpzAlg_check_install_runtime :: $(PppmpzAlg_check_install_runtime_dependencies)
PppmpzAlg_check_install_runtime ::
	$(silent) /cvmfs/bes3.ihep.ac.cn/bes3sw/Boss/6.6.5.p01/BesPolicy/BesPolicy-01-05-05/cmt/bes_check_installations.sh -files= -s=../share *.txt   -installdir=/workfs/bes/leizh/workarea-6.6.5.p01/InstallArea/share

endif

install ::
uninstall ::

#-- end of cmt_action_runner_header -----------------
#-- start of cleanup_header --------------

clean :: PppmpzAlg_check_install_runtimeclean ;
#	@cd .

ifndef PEDANTIC
.DEFAULT::
	$(echo) "(PppmpzAlg_check_install_runtime.make) $@: No rule for such target" >&2
else
.DEFAULT::
	$(error PEDANTIC: $@: No rule for such target)
endif

PppmpzAlg_check_install_runtimeclean ::
#-- end of cleanup_header ---------------
