#-- start of make_header -----------------

#====================================
#  Document omega_check_install_runtime
#
#   Generated Thu Aug  8 23:31:47 2019  by leizh
#
#====================================

include ${CMTROOT}/src/Makefile.core

ifdef tag
CMTEXTRATAGS = $(tag)
else
tag       = $(CMTCONFIG)
endif

cmt_omega_check_install_runtime_has_no_target_tag = 1

#--------------------------------------------------------

ifdef cmt_omega_check_install_runtime_has_target_tag

tags      = $(tag),$(CMTEXTRATAGS),target_omega_check_install_runtime

omega_tag = $(tag)

#cmt_local_tagfile_omega_check_install_runtime = $(omega_tag)_omega_check_install_runtime.make
cmt_local_tagfile_omega_check_install_runtime = $(bin)$(omega_tag)_omega_check_install_runtime.make

else

tags      = $(tag),$(CMTEXTRATAGS)

omega_tag = $(tag)

#cmt_local_tagfile_omega_check_install_runtime = $(omega_tag).make
cmt_local_tagfile_omega_check_install_runtime = $(bin)$(omega_tag).make

endif

include $(cmt_local_tagfile_omega_check_install_runtime)
#-include $(cmt_local_tagfile_omega_check_install_runtime)

ifdef cmt_omega_check_install_runtime_has_target_tag

cmt_final_setup_omega_check_install_runtime = $(bin)setup_omega_check_install_runtime.make
cmt_dependencies_in_omega_check_install_runtime = $(bin)dependencies_omega_check_install_runtime.in
#cmt_final_setup_omega_check_install_runtime = $(bin)omega_omega_check_install_runtimesetup.make
cmt_local_omega_check_install_runtime_makefile = $(bin)omega_check_install_runtime.make

else

cmt_final_setup_omega_check_install_runtime = $(bin)setup.make
cmt_dependencies_in_omega_check_install_runtime = $(bin)dependencies.in
#cmt_final_setup_omega_check_install_runtime = $(bin)omegasetup.make
cmt_local_omega_check_install_runtime_makefile = $(bin)omega_check_install_runtime.make

endif

#cmt_final_setup = $(bin)setup.make
#cmt_final_setup = $(bin)omegasetup.make

#omega_check_install_runtime :: ;

dirs ::
	@if test ! -r requirements ; then echo "No requirements file" ; fi; \
	  if test ! -d $(bin) ; then $(mkdir) -p $(bin) ; fi

javadirs ::
	@if test ! -d $(javabin) ; then $(mkdir) -p $(javabin) ; fi

srcdirs ::
	@if test ! -d $(src) ; then $(mkdir) -p $(src) ; fi

help ::
	$(echo) 'omega_check_install_runtime'

binobj = 
ifdef STRUCTURED_OUTPUT
binobj = omega_check_install_runtime/
#omega_check_install_runtime::
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
omega_check_install_runtime_once = 1
endif

ifdef omega_check_install_runtime_once

omega_check_install_runtimeactionstamp = $(bin)omega_check_install_runtime.actionstamp
#omega_check_install_runtimeactionstamp = omega_check_install_runtime.actionstamp

omega_check_install_runtime :: $(omega_check_install_runtimeactionstamp)
	$(echo) "omega_check_install_runtime ok"
#	@echo omega_check_install_runtime ok

#$(omega_check_install_runtimeactionstamp) :: $(omega_check_install_runtime_dependencies)
$(omega_check_install_runtimeactionstamp) ::
	$(silent) /afs/ihep.ac.cn/bes3/offline/Boss/6.6.5.p01/BesPolicy/BesPolicy-01-05-05/cmt/bes_check_installations.sh -files= -s=../share *.txt   -installdir=/share
	$(silent) cat /dev/null > $(omega_check_install_runtimeactionstamp)
#	@echo ok > $(omega_check_install_runtimeactionstamp)

omega_check_install_runtimeclean ::
	$(cleanup_silent) /bin/rm -f $(omega_check_install_runtimeactionstamp)

else

#omega_check_install_runtime :: $(omega_check_install_runtime_dependencies)
omega_check_install_runtime ::
	$(silent) /afs/ihep.ac.cn/bes3/offline/Boss/6.6.5.p01/BesPolicy/BesPolicy-01-05-05/cmt/bes_check_installations.sh -files= -s=../share *.txt   -installdir=/share

endif

install ::
uninstall ::

#-- end of cmt_action_runner_header -----------------
#-- start of cleanup_header --------------

clean :: omega_check_install_runtimeclean ;
#	@cd .

ifndef PEDANTIC
.DEFAULT::
	$(echo) "(omega_check_install_runtime.make) $@: No rule for such target" >&2
else
.DEFAULT::
	$(error PEDANTIC: $@: No rule for such target)
endif

omega_check_install_runtimeclean ::
#-- end of cleanup_header ---------------
