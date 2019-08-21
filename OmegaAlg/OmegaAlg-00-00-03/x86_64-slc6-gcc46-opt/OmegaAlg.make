#-- start of make_header -----------------

#====================================
#  Library OmegaAlg
#
#   Generated Wed Aug 21 17:14:29 2019  by leizh
#
#====================================

include ${CMTROOT}/src/Makefile.core

ifdef tag
CMTEXTRATAGS = $(tag)
else
tag       = $(CMTCONFIG)
endif

cmt_OmegaAlg_has_no_target_tag = 1

#--------------------------------------------------------

ifdef cmt_OmegaAlg_has_target_tag

tags      = $(tag),$(CMTEXTRATAGS),target_OmegaAlg

OmegaAlg_tag = $(tag)

#cmt_local_tagfile_OmegaAlg = $(OmegaAlg_tag)_OmegaAlg.make
cmt_local_tagfile_OmegaAlg = $(bin)$(OmegaAlg_tag)_OmegaAlg.make

else

tags      = $(tag),$(CMTEXTRATAGS)

OmegaAlg_tag = $(tag)

#cmt_local_tagfile_OmegaAlg = $(OmegaAlg_tag).make
cmt_local_tagfile_OmegaAlg = $(bin)$(OmegaAlg_tag).make

endif

include $(cmt_local_tagfile_OmegaAlg)
#-include $(cmt_local_tagfile_OmegaAlg)

ifdef cmt_OmegaAlg_has_target_tag

cmt_final_setup_OmegaAlg = $(bin)setup_OmegaAlg.make
cmt_dependencies_in_OmegaAlg = $(bin)dependencies_OmegaAlg.in
#cmt_final_setup_OmegaAlg = $(bin)OmegaAlg_OmegaAlgsetup.make
cmt_local_OmegaAlg_makefile = $(bin)OmegaAlg.make

else

cmt_final_setup_OmegaAlg = $(bin)setup.make
cmt_dependencies_in_OmegaAlg = $(bin)dependencies.in
#cmt_final_setup_OmegaAlg = $(bin)OmegaAlgsetup.make
cmt_local_OmegaAlg_makefile = $(bin)OmegaAlg.make

endif

#cmt_final_setup = $(bin)setup.make
#cmt_final_setup = $(bin)OmegaAlgsetup.make

#OmegaAlg :: ;

dirs ::
	@if test ! -r requirements ; then echo "No requirements file" ; fi; \
	  if test ! -d $(bin) ; then $(mkdir) -p $(bin) ; fi

javadirs ::
	@if test ! -d $(javabin) ; then $(mkdir) -p $(javabin) ; fi

srcdirs ::
	@if test ! -d $(src) ; then $(mkdir) -p $(src) ; fi

help ::
	$(echo) 'OmegaAlg'

binobj = 
ifdef STRUCTURED_OUTPUT
binobj = OmegaAlg/
#OmegaAlg::
#	@if test ! -d $(bin)$(binobj) ; then $(mkdir) -p $(bin)$(binobj) ; fi
#	$(echo) "STRUCTURED_OUTPUT="$(bin)$(binobj)
endif

${CMTROOT}/src/Makefile.core : ;
ifdef use_requirements
$(use_requirements) : ;
endif

#-- end of make_header ------------------
#-- start of libary_header ---------------

OmegaAlglibname   = $(bin)$(library_prefix)OmegaAlg$(library_suffix)
OmegaAlglib       = $(OmegaAlglibname).a
OmegaAlgstamp     = $(bin)OmegaAlg.stamp
OmegaAlgshstamp   = $(bin)OmegaAlg.shstamp

OmegaAlg :: dirs  OmegaAlgLIB
	$(echo) "OmegaAlg ok"

#-- end of libary_header ----------------

OmegaAlgLIB :: $(OmegaAlglib) $(OmegaAlgshstamp)
	@/bin/echo "------> OmegaAlg : library ok"

$(OmegaAlglib) :: $(bin)Omega.o $(bin)Omega_load.o $(bin)Omega_entries.o
	$(lib_echo) library
	$(lib_silent) cd $(bin); \
	  $(ar) $(OmegaAlglib) $?
	$(lib_silent) $(ranlib) $(OmegaAlglib)
	$(lib_silent) cat /dev/null >$(OmegaAlgstamp)

#------------------------------------------------------------------
#  Future improvement? to empty the object files after
#  storing in the library
#
##	  for f in $?; do \
##	    rm $${f}; touch $${f}; \
##	  done
#------------------------------------------------------------------

$(OmegaAlglibname).$(shlibsuffix) :: $(OmegaAlglib) $(OmegaAlgstamps)
	$(lib_silent) cd $(bin); QUIET=$(QUIET); $(make_shlib) "$(tags)" OmegaAlg $(OmegaAlg_shlibflags)

$(OmegaAlgshstamp) :: $(OmegaAlglibname).$(shlibsuffix)
	@if test -f $(OmegaAlglibname).$(shlibsuffix) ; then cat /dev/null >$(OmegaAlgshstamp) ; fi

OmegaAlgclean ::
	$(cleanup_echo) objects
	$(cleanup_silent) cd $(bin); /bin/rm -f $(bin)Omega.o $(bin)Omega_load.o $(bin)Omega_entries.o

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
OmegaAlginstallname = $(library_prefix)OmegaAlg$(library_suffix).$(shlibsuffix)

OmegaAlg :: OmegaAlginstall

install :: OmegaAlginstall

OmegaAlginstall :: $(install_dir)/$(OmegaAlginstallname)
	@if test ! "${installarea}" = ""; then\
	  echo "installation done"; \
	fi

$(install_dir)/$(OmegaAlginstallname) :: $(bin)$(OmegaAlginstallname)
	@if test ! "${installarea}" = ""; then \
	  cd $(bin); \
	  if test ! "$(install_dir)" = ""; then \
	    if test ! -d "$(install_dir)"; then \
	      mkdir -p $(install_dir); \
	    fi ; \
	    if test -d "$(install_dir)"; then \
	      echo "Installing library $(OmegaAlginstallname) into $(install_dir)"; \
	      if test -e $(install_dir)/$(OmegaAlginstallname); then \
	        $(cmt_uninstall_area_command) $(install_dir)/$(OmegaAlginstallname); \
	        $(cmt_uninstall_area_command) $(install_dir)/$(OmegaAlginstallname).cmtref; \
	      fi; \
	      $(cmt_install_area_command) `pwd`/$(OmegaAlginstallname) $(install_dir)/$(OmegaAlginstallname); \
	      echo `pwd`/$(OmegaAlginstallname) >$(install_dir)/$(OmegaAlginstallname).cmtref; \
	    fi \
          else \
	    echo "Cannot install library $(OmegaAlginstallname), no installation directory specified"; \
	  fi; \
	fi

OmegaAlgclean :: OmegaAlguninstall

uninstall :: OmegaAlguninstall

OmegaAlguninstall ::
	@if test ! "${installarea}" = ""; then \
	  cd $(bin); \
	  if test ! "$(install_dir)" = ""; then \
	    if test -d "$(install_dir)"; then \
	      echo "Removing installed library $(OmegaAlginstallname) from $(install_dir)"; \
	      $(cmt_uninstall_area_command) $(install_dir)/$(OmegaAlginstallname); \
	      $(cmt_uninstall_area_command) $(install_dir)/$(OmegaAlginstallname).cmtref; \
	    fi \
          else \
	    echo "Cannot uninstall library $(OmegaAlginstallname), no installation directory specified"; \
	  fi; \
	fi




#-- start of cpp_library -----------------

ifneq (-MMD -MP -MF $*.d -MQ $@,)

ifneq ($(MAKECMDGOALS),OmegaAlgclean)
ifneq ($(MAKECMDGOALS),uninstall)
-include $(bin)$(binobj)Omega.d

$(bin)$(binobj)Omega.d :

$(bin)$(binobj)Omega.o : $(cmt_final_setup_OmegaAlg)

$(bin)$(binobj)Omega.o : $(src)Omega.cxx
	$(cpp_echo) $(src)Omega.cxx
	$(cpp_silent) $(cppcomp) -MMD -MP -MF $*.d -MQ $@ -o $@ $(use_pp_cppflags) $(OmegaAlg_pp_cppflags) $(lib_OmegaAlg_pp_cppflags) $(Omega_pp_cppflags) $(use_cppflags) $(OmegaAlg_cppflags) $(lib_OmegaAlg_cppflags) $(Omega_cppflags) $(Omega_cxx_cppflags)  $(src)Omega.cxx
endif
endif

else
$(bin)OmegaAlg_dependencies.make : $(Omega_cxx_dependencies)

$(bin)OmegaAlg_dependencies.make : $(src)Omega.cxx

$(bin)$(binobj)Omega.o : $(Omega_cxx_dependencies)
	$(cpp_echo) $(src)Omega.cxx
	$(cpp_silent) $(cppcomp) -o $@ $(use_pp_cppflags) $(OmegaAlg_pp_cppflags) $(lib_OmegaAlg_pp_cppflags) $(Omega_pp_cppflags) $(use_cppflags) $(OmegaAlg_cppflags) $(lib_OmegaAlg_cppflags) $(Omega_cppflags) $(Omega_cxx_cppflags)  $(src)Omega.cxx

endif

#-- end of cpp_library ------------------
#-- start of cpp_library -----------------

ifneq (-MMD -MP -MF $*.d -MQ $@,)

ifneq ($(MAKECMDGOALS),OmegaAlgclean)
ifneq ($(MAKECMDGOALS),uninstall)
-include $(bin)$(binobj)Omega_load.d

$(bin)$(binobj)Omega_load.d :

$(bin)$(binobj)Omega_load.o : $(cmt_final_setup_OmegaAlg)

$(bin)$(binobj)Omega_load.o : $(src)components/Omega_load.cxx
	$(cpp_echo) $(src)components/Omega_load.cxx
	$(cpp_silent) $(cppcomp) -MMD -MP -MF $*.d -MQ $@ -o $@ $(use_pp_cppflags) $(OmegaAlg_pp_cppflags) $(lib_OmegaAlg_pp_cppflags) $(Omega_load_pp_cppflags) $(use_cppflags) $(OmegaAlg_cppflags) $(lib_OmegaAlg_cppflags) $(Omega_load_cppflags) $(Omega_load_cxx_cppflags) -I../src/components $(src)components/Omega_load.cxx
endif
endif

else
$(bin)OmegaAlg_dependencies.make : $(Omega_load_cxx_dependencies)

$(bin)OmegaAlg_dependencies.make : $(src)components/Omega_load.cxx

$(bin)$(binobj)Omega_load.o : $(Omega_load_cxx_dependencies)
	$(cpp_echo) $(src)components/Omega_load.cxx
	$(cpp_silent) $(cppcomp) -o $@ $(use_pp_cppflags) $(OmegaAlg_pp_cppflags) $(lib_OmegaAlg_pp_cppflags) $(Omega_load_pp_cppflags) $(use_cppflags) $(OmegaAlg_cppflags) $(lib_OmegaAlg_cppflags) $(Omega_load_cppflags) $(Omega_load_cxx_cppflags) -I../src/components $(src)components/Omega_load.cxx

endif

#-- end of cpp_library ------------------
#-- start of cpp_library -----------------

ifneq (-MMD -MP -MF $*.d -MQ $@,)

ifneq ($(MAKECMDGOALS),OmegaAlgclean)
ifneq ($(MAKECMDGOALS),uninstall)
-include $(bin)$(binobj)Omega_entries.d

$(bin)$(binobj)Omega_entries.d :

$(bin)$(binobj)Omega_entries.o : $(cmt_final_setup_OmegaAlg)

$(bin)$(binobj)Omega_entries.o : $(src)components/Omega_entries.cxx
	$(cpp_echo) $(src)components/Omega_entries.cxx
	$(cpp_silent) $(cppcomp) -MMD -MP -MF $*.d -MQ $@ -o $@ $(use_pp_cppflags) $(OmegaAlg_pp_cppflags) $(lib_OmegaAlg_pp_cppflags) $(Omega_entries_pp_cppflags) $(use_cppflags) $(OmegaAlg_cppflags) $(lib_OmegaAlg_cppflags) $(Omega_entries_cppflags) $(Omega_entries_cxx_cppflags) -I../src/components $(src)components/Omega_entries.cxx
endif
endif

else
$(bin)OmegaAlg_dependencies.make : $(Omega_entries_cxx_dependencies)

$(bin)OmegaAlg_dependencies.make : $(src)components/Omega_entries.cxx

$(bin)$(binobj)Omega_entries.o : $(Omega_entries_cxx_dependencies)
	$(cpp_echo) $(src)components/Omega_entries.cxx
	$(cpp_silent) $(cppcomp) -o $@ $(use_pp_cppflags) $(OmegaAlg_pp_cppflags) $(lib_OmegaAlg_pp_cppflags) $(Omega_entries_pp_cppflags) $(use_cppflags) $(OmegaAlg_cppflags) $(lib_OmegaAlg_cppflags) $(Omega_entries_cppflags) $(Omega_entries_cxx_cppflags) -I../src/components $(src)components/Omega_entries.cxx

endif

#-- end of cpp_library ------------------
#-- start of cleanup_header --------------

clean :: OmegaAlgclean ;
#	@cd .

ifndef PEDANTIC
.DEFAULT::
	$(echo) "(OmegaAlg.make) $@: No rule for such target" >&2
else
.DEFAULT::
	$(error PEDANTIC: $@: No rule for such target)
endif

OmegaAlgclean ::
#-- end of cleanup_header ---------------
#-- start of cleanup_library -------------
	$(cleanup_echo) library OmegaAlg
	-$(cleanup_silent) cd $(bin); /bin/rm -f $(library_prefix)OmegaAlg$(library_suffix).a $(library_prefix)OmegaAlg$(library_suffix).s? OmegaAlg.stamp OmegaAlg.shstamp
#-- end of cleanup_library ---------------
