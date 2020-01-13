#-- start of make_header -----------------

#====================================
#  Library PmAlg
#
#   Generated Tue Jan  7 09:21:11 2020  by leizh
#
#====================================

include ${CMTROOT}/src/Makefile.core

ifdef tag
CMTEXTRATAGS = $(tag)
else
tag       = $(CMTCONFIG)
endif

cmt_PmAlg_has_no_target_tag = 1

#--------------------------------------------------------

ifdef cmt_PmAlg_has_target_tag

tags      = $(tag),$(CMTEXTRATAGS),target_PmAlg

PmAlg_tag = $(tag)

#cmt_local_tagfile_PmAlg = $(PmAlg_tag)_PmAlg.make
cmt_local_tagfile_PmAlg = $(bin)$(PmAlg_tag)_PmAlg.make

else

tags      = $(tag),$(CMTEXTRATAGS)

PmAlg_tag = $(tag)

#cmt_local_tagfile_PmAlg = $(PmAlg_tag).make
cmt_local_tagfile_PmAlg = $(bin)$(PmAlg_tag).make

endif

include $(cmt_local_tagfile_PmAlg)
#-include $(cmt_local_tagfile_PmAlg)

ifdef cmt_PmAlg_has_target_tag

cmt_final_setup_PmAlg = $(bin)setup_PmAlg.make
cmt_dependencies_in_PmAlg = $(bin)dependencies_PmAlg.in
#cmt_final_setup_PmAlg = $(bin)PmAlg_PmAlgsetup.make
cmt_local_PmAlg_makefile = $(bin)PmAlg.make

else

cmt_final_setup_PmAlg = $(bin)setup.make
cmt_dependencies_in_PmAlg = $(bin)dependencies.in
#cmt_final_setup_PmAlg = $(bin)PmAlgsetup.make
cmt_local_PmAlg_makefile = $(bin)PmAlg.make

endif

#cmt_final_setup = $(bin)setup.make
#cmt_final_setup = $(bin)PmAlgsetup.make

#PmAlg :: ;

dirs ::
	@if test ! -r requirements ; then echo "No requirements file" ; fi; \
	  if test ! -d $(bin) ; then $(mkdir) -p $(bin) ; fi

javadirs ::
	@if test ! -d $(javabin) ; then $(mkdir) -p $(javabin) ; fi

srcdirs ::
	@if test ! -d $(src) ; then $(mkdir) -p $(src) ; fi

help ::
	$(echo) 'PmAlg'

binobj = 
ifdef STRUCTURED_OUTPUT
binobj = PmAlg/
#PmAlg::
#	@if test ! -d $(bin)$(binobj) ; then $(mkdir) -p $(bin)$(binobj) ; fi
#	$(echo) "STRUCTURED_OUTPUT="$(bin)$(binobj)
endif

${CMTROOT}/src/Makefile.core : ;
ifdef use_requirements
$(use_requirements) : ;
endif

#-- end of make_header ------------------
#-- start of libary_header ---------------

PmAlglibname   = $(bin)$(library_prefix)PmAlg$(library_suffix)
PmAlglib       = $(PmAlglibname).a
PmAlgstamp     = $(bin)PmAlg.stamp
PmAlgshstamp   = $(bin)PmAlg.shstamp

PmAlg :: dirs  PmAlgLIB
	$(echo) "PmAlg ok"

#-- end of libary_header ----------------

PmAlgLIB :: $(PmAlglib) $(PmAlgshstamp)
	@/bin/echo "------> PmAlg : library ok"

$(PmAlglib) :: $(bin)Pm.o $(bin)Pm_entries.o $(bin)Pm_load.o
	$(lib_echo) library
	$(lib_silent) cd $(bin); \
	  $(ar) $(PmAlglib) $?
	$(lib_silent) $(ranlib) $(PmAlglib)
	$(lib_silent) cat /dev/null >$(PmAlgstamp)

#------------------------------------------------------------------
#  Future improvement? to empty the object files after
#  storing in the library
#
##	  for f in $?; do \
##	    rm $${f}; touch $${f}; \
##	  done
#------------------------------------------------------------------

$(PmAlglibname).$(shlibsuffix) :: $(PmAlglib) $(PmAlgstamps)
	$(lib_silent) cd $(bin); QUIET=$(QUIET); $(make_shlib) "$(tags)" PmAlg $(PmAlg_shlibflags)

$(PmAlgshstamp) :: $(PmAlglibname).$(shlibsuffix)
	@if test -f $(PmAlglibname).$(shlibsuffix) ; then cat /dev/null >$(PmAlgshstamp) ; fi

PmAlgclean ::
	$(cleanup_echo) objects
	$(cleanup_silent) cd $(bin); /bin/rm -f $(bin)Pm.o $(bin)Pm_entries.o $(bin)Pm_load.o

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
PmAlginstallname = $(library_prefix)PmAlg$(library_suffix).$(shlibsuffix)

PmAlg :: PmAlginstall

install :: PmAlginstall

PmAlginstall :: $(install_dir)/$(PmAlginstallname)
	@if test ! "${installarea}" = ""; then\
	  echo "installation done"; \
	fi

$(install_dir)/$(PmAlginstallname) :: $(bin)$(PmAlginstallname)
	@if test ! "${installarea}" = ""; then \
	  cd $(bin); \
	  if test ! "$(install_dir)" = ""; then \
	    if test ! -d "$(install_dir)"; then \
	      mkdir -p $(install_dir); \
	    fi ; \
	    if test -d "$(install_dir)"; then \
	      echo "Installing library $(PmAlginstallname) into $(install_dir)"; \
	      if test -e $(install_dir)/$(PmAlginstallname); then \
	        $(cmt_uninstall_area_command) $(install_dir)/$(PmAlginstallname); \
	        $(cmt_uninstall_area_command) $(install_dir)/$(PmAlginstallname).cmtref; \
	      fi; \
	      $(cmt_install_area_command) `pwd`/$(PmAlginstallname) $(install_dir)/$(PmAlginstallname); \
	      echo `pwd`/$(PmAlginstallname) >$(install_dir)/$(PmAlginstallname).cmtref; \
	    fi \
          else \
	    echo "Cannot install library $(PmAlginstallname), no installation directory specified"; \
	  fi; \
	fi

PmAlgclean :: PmAlguninstall

uninstall :: PmAlguninstall

PmAlguninstall ::
	@if test ! "${installarea}" = ""; then \
	  cd $(bin); \
	  if test ! "$(install_dir)" = ""; then \
	    if test -d "$(install_dir)"; then \
	      echo "Removing installed library $(PmAlginstallname) from $(install_dir)"; \
	      $(cmt_uninstall_area_command) $(install_dir)/$(PmAlginstallname); \
	      $(cmt_uninstall_area_command) $(install_dir)/$(PmAlginstallname).cmtref; \
	    fi \
          else \
	    echo "Cannot uninstall library $(PmAlginstallname), no installation directory specified"; \
	  fi; \
	fi




#-- start of cpp_library -----------------

ifneq (-MMD -MP -MF $*.d -MQ $@,)

ifneq ($(MAKECMDGOALS),PmAlgclean)
ifneq ($(MAKECMDGOALS),uninstall)
-include $(bin)$(binobj)Pm.d

$(bin)$(binobj)Pm.d :

$(bin)$(binobj)Pm.o : $(cmt_final_setup_PmAlg)

$(bin)$(binobj)Pm.o : $(src)Pm.cxx
	$(cpp_echo) $(src)Pm.cxx
	$(cpp_silent) $(cppcomp) -MMD -MP -MF $*.d -MQ $@ -o $@ $(use_pp_cppflags) $(PmAlg_pp_cppflags) $(lib_PmAlg_pp_cppflags) $(Pm_pp_cppflags) $(use_cppflags) $(PmAlg_cppflags) $(lib_PmAlg_cppflags) $(Pm_cppflags) $(Pm_cxx_cppflags)  $(src)Pm.cxx
endif
endif

else
$(bin)PmAlg_dependencies.make : $(Pm_cxx_dependencies)

$(bin)PmAlg_dependencies.make : $(src)Pm.cxx

$(bin)$(binobj)Pm.o : $(Pm_cxx_dependencies)
	$(cpp_echo) $(src)Pm.cxx
	$(cpp_silent) $(cppcomp) -o $@ $(use_pp_cppflags) $(PmAlg_pp_cppflags) $(lib_PmAlg_pp_cppflags) $(Pm_pp_cppflags) $(use_cppflags) $(PmAlg_cppflags) $(lib_PmAlg_cppflags) $(Pm_cppflags) $(Pm_cxx_cppflags)  $(src)Pm.cxx

endif

#-- end of cpp_library ------------------
#-- start of cpp_library -----------------

ifneq (-MMD -MP -MF $*.d -MQ $@,)

ifneq ($(MAKECMDGOALS),PmAlgclean)
ifneq ($(MAKECMDGOALS),uninstall)
-include $(bin)$(binobj)Pm_entries.d

$(bin)$(binobj)Pm_entries.d :

$(bin)$(binobj)Pm_entries.o : $(cmt_final_setup_PmAlg)

$(bin)$(binobj)Pm_entries.o : $(src)components/Pm_entries.cxx
	$(cpp_echo) $(src)components/Pm_entries.cxx
	$(cpp_silent) $(cppcomp) -MMD -MP -MF $*.d -MQ $@ -o $@ $(use_pp_cppflags) $(PmAlg_pp_cppflags) $(lib_PmAlg_pp_cppflags) $(Pm_entries_pp_cppflags) $(use_cppflags) $(PmAlg_cppflags) $(lib_PmAlg_cppflags) $(Pm_entries_cppflags) $(Pm_entries_cxx_cppflags) -I../src/components $(src)components/Pm_entries.cxx
endif
endif

else
$(bin)PmAlg_dependencies.make : $(Pm_entries_cxx_dependencies)

$(bin)PmAlg_dependencies.make : $(src)components/Pm_entries.cxx

$(bin)$(binobj)Pm_entries.o : $(Pm_entries_cxx_dependencies)
	$(cpp_echo) $(src)components/Pm_entries.cxx
	$(cpp_silent) $(cppcomp) -o $@ $(use_pp_cppflags) $(PmAlg_pp_cppflags) $(lib_PmAlg_pp_cppflags) $(Pm_entries_pp_cppflags) $(use_cppflags) $(PmAlg_cppflags) $(lib_PmAlg_cppflags) $(Pm_entries_cppflags) $(Pm_entries_cxx_cppflags) -I../src/components $(src)components/Pm_entries.cxx

endif

#-- end of cpp_library ------------------
#-- start of cpp_library -----------------

ifneq (-MMD -MP -MF $*.d -MQ $@,)

ifneq ($(MAKECMDGOALS),PmAlgclean)
ifneq ($(MAKECMDGOALS),uninstall)
-include $(bin)$(binobj)Pm_load.d

$(bin)$(binobj)Pm_load.d :

$(bin)$(binobj)Pm_load.o : $(cmt_final_setup_PmAlg)

$(bin)$(binobj)Pm_load.o : $(src)components/Pm_load.cxx
	$(cpp_echo) $(src)components/Pm_load.cxx
	$(cpp_silent) $(cppcomp) -MMD -MP -MF $*.d -MQ $@ -o $@ $(use_pp_cppflags) $(PmAlg_pp_cppflags) $(lib_PmAlg_pp_cppflags) $(Pm_load_pp_cppflags) $(use_cppflags) $(PmAlg_cppflags) $(lib_PmAlg_cppflags) $(Pm_load_cppflags) $(Pm_load_cxx_cppflags) -I../src/components $(src)components/Pm_load.cxx
endif
endif

else
$(bin)PmAlg_dependencies.make : $(Pm_load_cxx_dependencies)

$(bin)PmAlg_dependencies.make : $(src)components/Pm_load.cxx

$(bin)$(binobj)Pm_load.o : $(Pm_load_cxx_dependencies)
	$(cpp_echo) $(src)components/Pm_load.cxx
	$(cpp_silent) $(cppcomp) -o $@ $(use_pp_cppflags) $(PmAlg_pp_cppflags) $(lib_PmAlg_pp_cppflags) $(Pm_load_pp_cppflags) $(use_cppflags) $(PmAlg_cppflags) $(lib_PmAlg_cppflags) $(Pm_load_cppflags) $(Pm_load_cxx_cppflags) -I../src/components $(src)components/Pm_load.cxx

endif

#-- end of cpp_library ------------------
#-- start of cleanup_header --------------

clean :: PmAlgclean ;
#	@cd .

ifndef PEDANTIC
.DEFAULT::
	$(echo) "(PmAlg.make) $@: No rule for such target" >&2
else
.DEFAULT::
	$(error PEDANTIC: $@: No rule for such target)
endif

PmAlgclean ::
#-- end of cleanup_header ---------------
#-- start of cleanup_library -------------
	$(cleanup_echo) library PmAlg
	-$(cleanup_silent) cd $(bin); /bin/rm -f $(library_prefix)PmAlg$(library_suffix).a $(library_prefix)PmAlg$(library_suffix).s? PmAlg.stamp PmAlg.shstamp
#-- end of cleanup_library ---------------
