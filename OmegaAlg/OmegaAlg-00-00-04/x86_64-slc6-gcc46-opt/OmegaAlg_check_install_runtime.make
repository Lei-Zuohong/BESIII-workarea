#-- start of make_header -----------------

#====================================
#  Document OmegaAlg_check_install_runtime
#
#   Generated Mon Aug 12 16:49:51 2019  by leizh
#
#====================================

include ${CMTROOT}/src/Makefile.core

ifdef tag
CMTEXTRATAGS = $(tag)
else
tag       = $(CMTCONFIG)
endif

cmt_OmegaAlg_check_install_runtime_has_no_target_tag = 1

#--------------------------------------------------------

ifdef cmt_OmegaAlg_check_install_runtime_has_target_tag

tags      = $(tag),$(CMTEXTRATAGS),target_OmegaAlg_check_install_runtime

OmegaAlg_tag = $(tag)

#cmt_local_tagfile_OmegaAlg_check_install_runtime = $(OmegaAlg_tag)_OmegaAlg_check_install_runtime.make
cmt_local_tagfile_OmegaAlg_check_install_runtime = $(bin)$(OmegaAlg_tag)_OmegaAlg_check_install_runtime.make

else

tags      = $(tag),$(CMTEXTRATAGS)

OmegaAlg_tag = $(tag)

#cmt_local_tagfile_OmegaAlg_check_install_runtime = $(OmegaAlg_tag).make
cmt_local_tagfile_OmegaAlg_check_install_runtime = $(bin)$(OmegaAlg_tag).make

endif

include $(cmt_local_tagfile_OmegaAlg_check_install_runtime)
#-include $(cmt_local_tagfile_OmegaAlg_check_install_runtime)

ifdef cmt_OmegaAlg_check_install_runtime_has_target_tag

cmt_final_setup_OmegaAlg_check_install_runtime = $(bin)setup_OmegaAlg_check_install_runtime.make
cmt_dependencies_in_OmegaAlg_check_install_runtime = $(bin)dependencies_OmegaAlg_check_install_runtime.in
#cmt_final_setup_OmegaAlg_check_install_runtime = $(bin)OmegaAlg_OmegaAlg_check_install_runtimesetup.make
cmt_local_OmegaAlg_check_install_runtime_makefile = $(bin)OmegaAlg_check_install_runtime.make

else

cmt_final_setup_OmegaAlg_check_install_runtime = $(bin)setup.make
cmt_dependencies_in_OmegaAlg_check_install_runtime = $(bin)dependencies.in
#cmt_final_setup_OmegaAlg_check_install_runtime = $(bin)OmegaAlgsetup.make
cmt_local_OmegaAlg_check_install_runtime_makefile = $(bin)OmegaAlg_check_install_runtime.make

endif

#cmt_final_setup = $(bin)setup.make
#cmt_final_setup = $(bin)OmegaAlgsetup.make

#OmegaAlg_check_install_runtime :: ;

dirs ::
	@if test ! -r requirements ; then echo "No requirements file" ; fi; \
	  if test ! -d $(bin) ; then $(mkdir) -p $(bin) ; fi

javadirs ::
	@if test ! -d $(javabin) ; then $(mkdir) -p $(javabin) ; fi

srcdirs ::
	@if test ! -d $(src) ; then $(mkdir) -p $(src) ; fi

help ::
	$(echo) 'OmegaAlg_check_install_runtime'

binobj = 
ifdef STRUCTURED_OUTPUT
binobj = OmegaAlg_check_install_runtime/
#OmegaAlg_check_install_runtime::
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
OmegaAlg_check_install_runtime_once = 1
endif

ifdef OmegaAlg_check_install_runtime_once

OmegaAlg_check_install_runtimeactionstamp = $(bin)OmegaAlg_check_install_runtime.actionstamp
#OmegaAlg_check_install_runtimeactionstamp = OmegaAlg_check_install_runtime.actionstamp

OmegaAlg_check_install_runtime :: $(OmegaAlg_check_install_runtimeactionstamp)
	$(echo) "OmegaAlg_check_install_runtime ok"
#	@echo OmegaAlg_check_install_runtime ok

#$(OmegaAlg_check_install_runtimeactionstamp) :: $(OmegaAlg_check_install_runtime_dependencies)
$(OmegaAlg_check_install_runtimeactionstamp) ::
	$(silent) /afs/ihep.ac.cn/bes3/offline/Boss/6.6.5.p01/BesPolicy/BesPolicy-01-05-05/cmt/bes_check_installations.sh -files= -s=../share *.txt   -installdir=/workfs/bes/leizh/workarea-6.6.5.p01/InstallArea/share
	$(silent) cat /dev/null > $(OmegaAlg_check_install_runtimeactionstamp)
#	@echo ok > $(OmegaAlg_check_install_runtimeactionstamp)

OmegaAlg_check_install_runtimeclean ::
	$(cleanup_silent) /bin/rm -f $(OmegaAlg_check_install_runtimeactionstamp)

else

#OmegaAlg_check_install_runtime :: $(OmegaAlg_check_install_runtime_dependencies)
OmegaAlg_check_install_runtime ::
	$(silent) /afs/ihep.ac.cn/bes3/offline/Boss/6.6.5.p01/BesPolicy/BesPolicy-01-05-05/cmt/bes_check_installations.sh -files= -s=../share *.txt   -installdir=/workfs/bes/leizh/workarea-6.6.5.p01/InstallArea/share

endif

install ::
uninstall ::

#-- end of cmt_action_runner_header -----------------
#-- start of cleanup_header --------------

clean :: OmegaAlg_check_install_runtimeclean ;
#	@cd .

ifndef PEDANTIC
.DEFAULT::
	$(echo) "(OmegaAlg_check_install_runtime.make) $@: No rule for such target" >&2
else
.DEFAULT::
	$(error PEDANTIC: $@: No rule for such target)
endif

OmegaAlg_check_install_runtimeclean ::
#-- end of cleanup_header ---------------
