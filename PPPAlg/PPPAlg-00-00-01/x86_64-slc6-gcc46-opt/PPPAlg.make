#-- start of make_header -----------------

#====================================
#  Library PPPAlg
#
#   Generated Wed Jul 15 14:22:20 2020  by leizh
#
#====================================

include ${CMTROOT}/src/Makefile.core

ifdef tag
CMTEXTRATAGS = $(tag)
else
tag       = $(CMTCONFIG)
endif

cmt_PPPAlg_has_no_target_tag = 1

#--------------------------------------------------------

ifdef cmt_PPPAlg_has_target_tag

tags      = $(tag),$(CMTEXTRATAGS),target_PPPAlg

PPPAlg_tag = $(tag)

#cmt_local_tagfile_PPPAlg = $(PPPAlg_tag)_PPPAlg.make
cmt_local_tagfile_PPPAlg = $(bin)$(PPPAlg_tag)_PPPAlg.make

else

tags      = $(tag),$(CMTEXTRATAGS)

PPPAlg_tag = $(tag)

#cmt_local_tagfile_PPPAlg = $(PPPAlg_tag).make
cmt_local_tagfile_PPPAlg = $(bin)$(PPPAlg_tag).make

endif

include $(cmt_local_tagfile_PPPAlg)
#-include $(cmt_local_tagfile_PPPAlg)

ifdef cmt_PPPAlg_has_target_tag

cmt_final_setup_PPPAlg = $(bin)setup_PPPAlg.make
cmt_dependencies_in_PPPAlg = $(bin)dependencies_PPPAlg.in
#cmt_final_setup_PPPAlg = $(bin)PPPAlg_PPPAlgsetup.make
cmt_local_PPPAlg_makefile = $(bin)PPPAlg.make

else

cmt_final_setup_PPPAlg = $(bin)setup.make
cmt_dependencies_in_PPPAlg = $(bin)dependencies.in
#cmt_final_setup_PPPAlg = $(bin)PPPAlgsetup.make
cmt_local_PPPAlg_makefile = $(bin)PPPAlg.make

endif

#cmt_final_setup = $(bin)setup.make
#cmt_final_setup = $(bin)PPPAlgsetup.make

#PPPAlg :: ;

dirs ::
	@if test ! -r requirements ; then echo "No requirements file" ; fi; \
	  if test ! -d $(bin) ; then $(mkdir) -p $(bin) ; fi

javadirs ::
	@if test ! -d $(javabin) ; then $(mkdir) -p $(javabin) ; fi

srcdirs ::
	@if test ! -d $(src) ; then $(mkdir) -p $(src) ; fi

help ::
	$(echo) 'PPPAlg'

binobj = 
ifdef STRUCTURED_OUTPUT
binobj = PPPAlg/
#PPPAlg::
#	@if test ! -d $(bin)$(binobj) ; then $(mkdir) -p $(bin)$(binobj) ; fi
#	$(echo) "STRUCTURED_OUTPUT="$(bin)$(binobj)
endif

${CMTROOT}/src/Makefile.core : ;
ifdef use_requirements
$(use_requirements) : ;
endif

#-- end of make_header ------------------
#-- start of libary_header ---------------

PPPAlglibname   = $(bin)$(library_prefix)PPPAlg$(library_suffix)
PPPAlglib       = $(PPPAlglibname).a
PPPAlgstamp     = $(bin)PPPAlg.stamp
PPPAlgshstamp   = $(bin)PPPAlg.shstamp

PPPAlg :: dirs  PPPAlgLIB
	$(echo) "PPPAlg ok"

#-- end of libary_header ----------------

PPPAlgLIB :: $(PPPAlglib) $(PPPAlgshstamp)
	@/bin/echo "------> PPPAlg : library ok"

$(PPPAlglib) :: $(bin)PPP.o $(bin)PPP_entries.o $(bin)PPP_load.o
	$(lib_echo) library
	$(lib_silent) cd $(bin); \
	  $(ar) $(PPPAlglib) $?
	$(lib_silent) $(ranlib) $(PPPAlglib)
	$(lib_silent) cat /dev/null >$(PPPAlgstamp)

#------------------------------------------------------------------
#  Future improvement? to empty the object files after
#  storing in the library
#
##	  for f in $?; do \
##	    rm $${f}; touch $${f}; \
##	  done
#------------------------------------------------------------------

$(PPPAlglibname).$(shlibsuffix) :: $(PPPAlglib) $(PPPAlgstamps)
	$(lib_silent) cd $(bin); QUIET=$(QUIET); $(make_shlib) "$(tags)" PPPAlg $(PPPAlg_shlibflags)

$(PPPAlgshstamp) :: $(PPPAlglibname).$(shlibsuffix)
	@if test -f $(PPPAlglibname).$(shlibsuffix) ; then cat /dev/null >$(PPPAlgshstamp) ; fi

PPPAlgclean ::
	$(cleanup_echo) objects
	$(cleanup_silent) cd $(bin); /bin/rm -f $(bin)PPP.o $(bin)PPP_entries.o $(bin)PPP_load.o

#-----------------------------------------------------------------
#
#  New section for automatic installation
#
#-----------------------------------------------------------------

ifeq ($(INSTALLAREA),)
installarea = $(CMTINSTALLAREA)
else
ifeq ($(findstring `,$(INSTALLAREA)),`)
installarea = $(shell $(subst `,, $(INSTALLAREA)))
else
installarea = $(INSTALLAREA)
endif
endif

install_dir = ${installarea}/${CMTCONFIG}/lib
PPPAlginstallname = $(library_prefix)PPPAlg$(library_suffix).$(shlibsuffix)

PPPAlg :: PPPAlginstall

install :: PPPAlginstall

PPPAlginstall :: $(install_dir)/$(PPPAlginstallname)
	@if test ! "${installarea}" = ""; then\
	  echo "installation done"; \
	fi

$(install_dir)/$(PPPAlginstallname) :: $(bin)$(PPPAlginstallname)
	@if test ! "${installarea}" = ""; then \
	  cd $(bin); \
	  if test ! "$(install_dir)" = ""; then \
	    if test ! -d "$(install_dir)"; then \
	      mkdir -p $(install_dir); \
	    fi ; \
	    if test -d "$(install_dir)"; then \
	      echo "Installing library $(PPPAlginstallname) into $(install_dir)"; \
	      if test -e $(install_dir)/$(PPPAlginstallname); then \
	        $(cmt_uninstall_area_command) $(install_dir)/$(PPPAlginstallname); \
	        $(cmt_uninstall_area_command) $(install_dir)/$(PPPAlginstallname).cmtref; \
	      fi; \
	      $(cmt_install_area_command) `pwd`/$(PPPAlginstallname) $(install_dir)/$(PPPAlginstallname); \
	      echo `pwd`/$(PPPAlginstallname) >$(install_dir)/$(PPPAlginstallname).cmtref; \
	    fi \
          else \
	    echo "Cannot install library $(PPPAlginstallname), no installation directory specified"; \
	  fi; \
	fi

PPPAlgclean :: PPPAlguninstall

uninstall :: PPPAlguninstall

PPPAlguninstall ::
	@if test ! "${installarea}" = ""; then \
	  cd $(bin); \
	  if test ! "$(install_dir)" = ""; then \
	    if test -d "$(install_dir)"; then \
	      echo "Removing installed library $(PPPAlginstallname) from $(install_dir)"; \
	      $(cmt_uninstall_area_command) $(install_dir)/$(PPPAlginstallname); \
	      $(cmt_uninstall_area_command) $(install_dir)/$(PPPAlginstallname).cmtref; \
	    fi \
          else \
	    echo "Cannot uninstall library $(PPPAlginstallname), no installation directory specified"; \
	  fi; \
	fi




#-- start of cpp_library -----------------

ifneq (-MMD -MP -MF $*.d -MQ $@,)

ifneq ($(MAKECMDGOALS),PPPAlgclean)
ifneq ($(MAKECMDGOALS),uninstall)
-include $(bin)$(binobj)PPP.d

$(bin)$(binobj)PPP.d :

$(bin)$(binobj)PPP.o : $(cmt_final_setup_PPPAlg)

$(bin)$(binobj)PPP.o : $(src)PPP.cxx
	$(cpp_echo) $(src)PPP.cxx
	$(cpp_silent) $(cppcomp) -MMD -MP -MF $*.d -MQ $@ -o $@ $(use_pp_cppflags) $(PPPAlg_pp_cppflags) $(lib_PPPAlg_pp_cppflags) $(PPP_pp_cppflags) $(use_cppflags) $(PPPAlg_cppflags) $(lib_PPPAlg_cppflags) $(PPP_cppflags) $(PPP_cxx_cppflags)  $(src)PPP.cxx
endif
endif

else
$(bin)PPPAlg_dependencies.make : $(PPP_cxx_dependencies)

$(bin)PPPAlg_dependencies.make : $(src)PPP.cxx

$(bin)$(binobj)PPP.o : $(PPP_cxx_dependencies)
	$(cpp_echo) $(src)PPP.cxx
	$(cpp_silent) $(cppcomp) -o $@ $(use_pp_cppflags) $(PPPAlg_pp_cppflags) $(lib_PPPAlg_pp_cppflags) $(PPP_pp_cppflags) $(use_cppflags) $(PPPAlg_cppflags) $(lib_PPPAlg_cppflags) $(PPP_cppflags) $(PPP_cxx_cppflags)  $(src)PPP.cxx

endif

#-- end of cpp_library ------------------
#-- start of cpp_library -----------------

ifneq (-MMD -MP -MF $*.d -MQ $@,)

ifneq ($(MAKECMDGOALS),PPPAlgclean)
ifneq ($(MAKECMDGOALS),uninstall)
-include $(bin)$(binobj)PPP_entries.d

$(bin)$(binobj)PPP_entries.d :

$(bin)$(binobj)PPP_entries.o : $(cmt_final_setup_PPPAlg)

$(bin)$(binobj)PPP_entries.o : $(src)components/PPP_entries.cxx
	$(cpp_echo) $(src)components/PPP_entries.cxx
	$(cpp_silent) $(cppcomp) -MMD -MP -MF $*.d -MQ $@ -o $@ $(use_pp_cppflags) $(PPPAlg_pp_cppflags) $(lib_PPPAlg_pp_cppflags) $(PPP_entries_pp_cppflags) $(use_cppflags) $(PPPAlg_cppflags) $(lib_PPPAlg_cppflags) $(PPP_entries_cppflags) $(PPP_entries_cxx_cppflags) -I../src/components $(src)components/PPP_entries.cxx
endif
endif

else
$(bin)PPPAlg_dependencies.make : $(PPP_entries_cxx_dependencies)

$(bin)PPPAlg_dependencies.make : $(src)components/PPP_entries.cxx

$(bin)$(binobj)PPP_entries.o : $(PPP_entries_cxx_dependencies)
	$(cpp_echo) $(src)components/PPP_entries.cxx
	$(cpp_silent) $(cppcomp) -o $@ $(use_pp_cppflags) $(PPPAlg_pp_cppflags) $(lib_PPPAlg_pp_cppflags) $(PPP_entries_pp_cppflags) $(use_cppflags) $(PPPAlg_cppflags) $(lib_PPPAlg_cppflags) $(PPP_entries_cppflags) $(PPP_entries_cxx_cppflags) -I../src/components $(src)components/PPP_entries.cxx

endif

#-- end of cpp_library ------------------
#-- start of cpp_library -----------------

ifneq (-MMD -MP -MF $*.d -MQ $@,)

ifneq ($(MAKECMDGOALS),PPPAlgclean)
ifneq ($(MAKECMDGOALS),uninstall)
-include $(bin)$(binobj)PPP_load.d

$(bin)$(binobj)PPP_load.d :

$(bin)$(binobj)PPP_load.o : $(cmt_final_setup_PPPAlg)

$(bin)$(binobj)PPP_load.o : $(src)components/PPP_load.cxx
	$(cpp_echo) $(src)components/PPP_load.cxx
	$(cpp_silent) $(cppcomp) -MMD -MP -MF $*.d -MQ $@ -o $@ $(use_pp_cppflags) $(PPPAlg_pp_cppflags) $(lib_PPPAlg_pp_cppflags) $(PPP_load_pp_cppflags) $(use_cppflags) $(PPPAlg_cppflags) $(lib_PPPAlg_cppflags) $(PPP_load_cppflags) $(PPP_load_cxx_cppflags) -I../src/components $(src)components/PPP_load.cxx
endif
endif

else
$(bin)PPPAlg_dependencies.make : $(PPP_load_cxx_dependencies)

$(bin)PPPAlg_dependencies.make : $(src)components/PPP_load.cxx

$(bin)$(binobj)PPP_load.o : $(PPP_load_cxx_dependencies)
	$(cpp_echo) $(src)components/PPP_load.cxx
	$(cpp_silent) $(cppcomp) -o $@ $(use_pp_cppflags) $(PPPAlg_pp_cppflags) $(lib_PPPAlg_pp_cppflags) $(PPP_load_pp_cppflags) $(use_cppflags) $(PPPAlg_cppflags) $(lib_PPPAlg_cppflags) $(PPP_load_cppflags) $(PPP_load_cxx_cppflags) -I../src/components $(src)components/PPP_load.cxx

endif

#-- end of cpp_library ------------------
#-- start of cleanup_header --------------

clean :: PPPAlgclean ;
#	@cd .

ifndef PEDANTIC
.DEFAULT::
	$(echo) "(PPPAlg.make) $@: No rule for such target" >&2
else
.DEFAULT::
	$(error PEDANTIC: $@: No rule for such target)
endif

PPPAlgclean ::
#-- end of cleanup_header ---------------
#-- start of cleanup_library -------------
	$(cleanup_echo) library PPPAlg
	-$(cleanup_silent) cd $(bin); /bin/rm -f $(library_prefix)PPPAlg$(library_suffix).a $(library_prefix)PPPAlg$(library_suffix).s? PPPAlg.stamp PPPAlg.shstamp
#-- end of cleanup_library ---------------
