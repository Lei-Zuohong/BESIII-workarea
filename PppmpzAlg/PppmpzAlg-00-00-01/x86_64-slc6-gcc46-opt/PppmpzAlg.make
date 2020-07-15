#-- start of make_header -----------------

#====================================
#  Library PppmpzAlg
#
#   Generated Wed Jul 15 15:26:02 2020  by leizh
#
#====================================

include ${CMTROOT}/src/Makefile.core

ifdef tag
CMTEXTRATAGS = $(tag)
else
tag       = $(CMTCONFIG)
endif

cmt_PppmpzAlg_has_no_target_tag = 1

#--------------------------------------------------------

ifdef cmt_PppmpzAlg_has_target_tag

tags      = $(tag),$(CMTEXTRATAGS),target_PppmpzAlg

PppmpzAlg_tag = $(tag)

#cmt_local_tagfile_PppmpzAlg = $(PppmpzAlg_tag)_PppmpzAlg.make
cmt_local_tagfile_PppmpzAlg = $(bin)$(PppmpzAlg_tag)_PppmpzAlg.make

else

tags      = $(tag),$(CMTEXTRATAGS)

PppmpzAlg_tag = $(tag)

#cmt_local_tagfile_PppmpzAlg = $(PppmpzAlg_tag).make
cmt_local_tagfile_PppmpzAlg = $(bin)$(PppmpzAlg_tag).make

endif

include $(cmt_local_tagfile_PppmpzAlg)
#-include $(cmt_local_tagfile_PppmpzAlg)

ifdef cmt_PppmpzAlg_has_target_tag

cmt_final_setup_PppmpzAlg = $(bin)setup_PppmpzAlg.make
cmt_dependencies_in_PppmpzAlg = $(bin)dependencies_PppmpzAlg.in
#cmt_final_setup_PppmpzAlg = $(bin)PppmpzAlg_PppmpzAlgsetup.make
cmt_local_PppmpzAlg_makefile = $(bin)PppmpzAlg.make

else

cmt_final_setup_PppmpzAlg = $(bin)setup.make
cmt_dependencies_in_PppmpzAlg = $(bin)dependencies.in
#cmt_final_setup_PppmpzAlg = $(bin)PppmpzAlgsetup.make
cmt_local_PppmpzAlg_makefile = $(bin)PppmpzAlg.make

endif

#cmt_final_setup = $(bin)setup.make
#cmt_final_setup = $(bin)PppmpzAlgsetup.make

#PppmpzAlg :: ;

dirs ::
	@if test ! -r requirements ; then echo "No requirements file" ; fi; \
	  if test ! -d $(bin) ; then $(mkdir) -p $(bin) ; fi

javadirs ::
	@if test ! -d $(javabin) ; then $(mkdir) -p $(javabin) ; fi

srcdirs ::
	@if test ! -d $(src) ; then $(mkdir) -p $(src) ; fi

help ::
	$(echo) 'PppmpzAlg'

binobj = 
ifdef STRUCTURED_OUTPUT
binobj = PppmpzAlg/
#PppmpzAlg::
#	@if test ! -d $(bin)$(binobj) ; then $(mkdir) -p $(bin)$(binobj) ; fi
#	$(echo) "STRUCTURED_OUTPUT="$(bin)$(binobj)
endif

${CMTROOT}/src/Makefile.core : ;
ifdef use_requirements
$(use_requirements) : ;
endif

#-- end of make_header ------------------
#-- start of libary_header ---------------

PppmpzAlglibname   = $(bin)$(library_prefix)PppmpzAlg$(library_suffix)
PppmpzAlglib       = $(PppmpzAlglibname).a
PppmpzAlgstamp     = $(bin)PppmpzAlg.stamp
PppmpzAlgshstamp   = $(bin)PppmpzAlg.shstamp

PppmpzAlg :: dirs  PppmpzAlgLIB
	$(echo) "PppmpzAlg ok"

#-- end of libary_header ----------------

PppmpzAlgLIB :: $(PppmpzAlglib) $(PppmpzAlgshstamp)
	@/bin/echo "------> PppmpzAlg : library ok"

$(PppmpzAlglib) :: $(bin)Pppmpz.o $(bin)Pppmpz_entries.o $(bin)Pppmpz_load.o
	$(lib_echo) library
	$(lib_silent) cd $(bin); \
	  $(ar) $(PppmpzAlglib) $?
	$(lib_silent) $(ranlib) $(PppmpzAlglib)
	$(lib_silent) cat /dev/null >$(PppmpzAlgstamp)

#------------------------------------------------------------------
#  Future improvement? to empty the object files after
#  storing in the library
#
##	  for f in $?; do \
##	    rm $${f}; touch $${f}; \
##	  done
#------------------------------------------------------------------

$(PppmpzAlglibname).$(shlibsuffix) :: $(PppmpzAlglib) $(PppmpzAlgstamps)
	$(lib_silent) cd $(bin); QUIET=$(QUIET); $(make_shlib) "$(tags)" PppmpzAlg $(PppmpzAlg_shlibflags)

$(PppmpzAlgshstamp) :: $(PppmpzAlglibname).$(shlibsuffix)
	@if test -f $(PppmpzAlglibname).$(shlibsuffix) ; then cat /dev/null >$(PppmpzAlgshstamp) ; fi

PppmpzAlgclean ::
	$(cleanup_echo) objects
	$(cleanup_silent) cd $(bin); /bin/rm -f $(bin)Pppmpz.o $(bin)Pppmpz_entries.o $(bin)Pppmpz_load.o

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
PppmpzAlginstallname = $(library_prefix)PppmpzAlg$(library_suffix).$(shlibsuffix)

PppmpzAlg :: PppmpzAlginstall

install :: PppmpzAlginstall

PppmpzAlginstall :: $(install_dir)/$(PppmpzAlginstallname)
	@if test ! "${installarea}" = ""; then\
	  echo "installation done"; \
	fi

$(install_dir)/$(PppmpzAlginstallname) :: $(bin)$(PppmpzAlginstallname)
	@if test ! "${installarea}" = ""; then \
	  cd $(bin); \
	  if test ! "$(install_dir)" = ""; then \
	    if test ! -d "$(install_dir)"; then \
	      mkdir -p $(install_dir); \
	    fi ; \
	    if test -d "$(install_dir)"; then \
	      echo "Installing library $(PppmpzAlginstallname) into $(install_dir)"; \
	      if test -e $(install_dir)/$(PppmpzAlginstallname); then \
	        $(cmt_uninstall_area_command) $(install_dir)/$(PppmpzAlginstallname); \
	        $(cmt_uninstall_area_command) $(install_dir)/$(PppmpzAlginstallname).cmtref; \
	      fi; \
	      $(cmt_install_area_command) `pwd`/$(PppmpzAlginstallname) $(install_dir)/$(PppmpzAlginstallname); \
	      echo `pwd`/$(PppmpzAlginstallname) >$(install_dir)/$(PppmpzAlginstallname).cmtref; \
	    fi \
          else \
	    echo "Cannot install library $(PppmpzAlginstallname), no installation directory specified"; \
	  fi; \
	fi

PppmpzAlgclean :: PppmpzAlguninstall

uninstall :: PppmpzAlguninstall

PppmpzAlguninstall ::
	@if test ! "${installarea}" = ""; then \
	  cd $(bin); \
	  if test ! "$(install_dir)" = ""; then \
	    if test -d "$(install_dir)"; then \
	      echo "Removing installed library $(PppmpzAlginstallname) from $(install_dir)"; \
	      $(cmt_uninstall_area_command) $(install_dir)/$(PppmpzAlginstallname); \
	      $(cmt_uninstall_area_command) $(install_dir)/$(PppmpzAlginstallname).cmtref; \
	    fi \
          else \
	    echo "Cannot uninstall library $(PppmpzAlginstallname), no installation directory specified"; \
	  fi; \
	fi




#-- start of cpp_library -----------------

ifneq (-MMD -MP -MF $*.d -MQ $@,)

ifneq ($(MAKECMDGOALS),PppmpzAlgclean)
ifneq ($(MAKECMDGOALS),uninstall)
-include $(bin)$(binobj)Pppmpz.d

$(bin)$(binobj)Pppmpz.d :

$(bin)$(binobj)Pppmpz.o : $(cmt_final_setup_PppmpzAlg)

$(bin)$(binobj)Pppmpz.o : $(src)Pppmpz.cxx
	$(cpp_echo) $(src)Pppmpz.cxx
	$(cpp_silent) $(cppcomp) -MMD -MP -MF $*.d -MQ $@ -o $@ $(use_pp_cppflags) $(PppmpzAlg_pp_cppflags) $(lib_PppmpzAlg_pp_cppflags) $(Pppmpz_pp_cppflags) $(use_cppflags) $(PppmpzAlg_cppflags) $(lib_PppmpzAlg_cppflags) $(Pppmpz_cppflags) $(Pppmpz_cxx_cppflags)  $(src)Pppmpz.cxx
endif
endif

else
$(bin)PppmpzAlg_dependencies.make : $(Pppmpz_cxx_dependencies)

$(bin)PppmpzAlg_dependencies.make : $(src)Pppmpz.cxx

$(bin)$(binobj)Pppmpz.o : $(Pppmpz_cxx_dependencies)
	$(cpp_echo) $(src)Pppmpz.cxx
	$(cpp_silent) $(cppcomp) -o $@ $(use_pp_cppflags) $(PppmpzAlg_pp_cppflags) $(lib_PppmpzAlg_pp_cppflags) $(Pppmpz_pp_cppflags) $(use_cppflags) $(PppmpzAlg_cppflags) $(lib_PppmpzAlg_cppflags) $(Pppmpz_cppflags) $(Pppmpz_cxx_cppflags)  $(src)Pppmpz.cxx

endif

#-- end of cpp_library ------------------
#-- start of cpp_library -----------------

ifneq (-MMD -MP -MF $*.d -MQ $@,)

ifneq ($(MAKECMDGOALS),PppmpzAlgclean)
ifneq ($(MAKECMDGOALS),uninstall)
-include $(bin)$(binobj)Pppmpz_entries.d

$(bin)$(binobj)Pppmpz_entries.d :

$(bin)$(binobj)Pppmpz_entries.o : $(cmt_final_setup_PppmpzAlg)

$(bin)$(binobj)Pppmpz_entries.o : $(src)components/Pppmpz_entries.cxx
	$(cpp_echo) $(src)components/Pppmpz_entries.cxx
	$(cpp_silent) $(cppcomp) -MMD -MP -MF $*.d -MQ $@ -o $@ $(use_pp_cppflags) $(PppmpzAlg_pp_cppflags) $(lib_PppmpzAlg_pp_cppflags) $(Pppmpz_entries_pp_cppflags) $(use_cppflags) $(PppmpzAlg_cppflags) $(lib_PppmpzAlg_cppflags) $(Pppmpz_entries_cppflags) $(Pppmpz_entries_cxx_cppflags) -I../src/components $(src)components/Pppmpz_entries.cxx
endif
endif

else
$(bin)PppmpzAlg_dependencies.make : $(Pppmpz_entries_cxx_dependencies)

$(bin)PppmpzAlg_dependencies.make : $(src)components/Pppmpz_entries.cxx

$(bin)$(binobj)Pppmpz_entries.o : $(Pppmpz_entries_cxx_dependencies)
	$(cpp_echo) $(src)components/Pppmpz_entries.cxx
	$(cpp_silent) $(cppcomp) -o $@ $(use_pp_cppflags) $(PppmpzAlg_pp_cppflags) $(lib_PppmpzAlg_pp_cppflags) $(Pppmpz_entries_pp_cppflags) $(use_cppflags) $(PppmpzAlg_cppflags) $(lib_PppmpzAlg_cppflags) $(Pppmpz_entries_cppflags) $(Pppmpz_entries_cxx_cppflags) -I../src/components $(src)components/Pppmpz_entries.cxx

endif

#-- end of cpp_library ------------------
#-- start of cpp_library -----------------

ifneq (-MMD -MP -MF $*.d -MQ $@,)

ifneq ($(MAKECMDGOALS),PppmpzAlgclean)
ifneq ($(MAKECMDGOALS),uninstall)
-include $(bin)$(binobj)Pppmpz_load.d

$(bin)$(binobj)Pppmpz_load.d :

$(bin)$(binobj)Pppmpz_load.o : $(cmt_final_setup_PppmpzAlg)

$(bin)$(binobj)Pppmpz_load.o : $(src)components/Pppmpz_load.cxx
	$(cpp_echo) $(src)components/Pppmpz_load.cxx
	$(cpp_silent) $(cppcomp) -MMD -MP -MF $*.d -MQ $@ -o $@ $(use_pp_cppflags) $(PppmpzAlg_pp_cppflags) $(lib_PppmpzAlg_pp_cppflags) $(Pppmpz_load_pp_cppflags) $(use_cppflags) $(PppmpzAlg_cppflags) $(lib_PppmpzAlg_cppflags) $(Pppmpz_load_cppflags) $(Pppmpz_load_cxx_cppflags) -I../src/components $(src)components/Pppmpz_load.cxx
endif
endif

else
$(bin)PppmpzAlg_dependencies.make : $(Pppmpz_load_cxx_dependencies)

$(bin)PppmpzAlg_dependencies.make : $(src)components/Pppmpz_load.cxx

$(bin)$(binobj)Pppmpz_load.o : $(Pppmpz_load_cxx_dependencies)
	$(cpp_echo) $(src)components/Pppmpz_load.cxx
	$(cpp_silent) $(cppcomp) -o $@ $(use_pp_cppflags) $(PppmpzAlg_pp_cppflags) $(lib_PppmpzAlg_pp_cppflags) $(Pppmpz_load_pp_cppflags) $(use_cppflags) $(PppmpzAlg_cppflags) $(lib_PppmpzAlg_cppflags) $(Pppmpz_load_cppflags) $(Pppmpz_load_cxx_cppflags) -I../src/components $(src)components/Pppmpz_load.cxx

endif

#-- end of cpp_library ------------------
#-- start of cleanup_header --------------

clean :: PppmpzAlgclean ;
#	@cd .

ifndef PEDANTIC
.DEFAULT::
	$(echo) "(PppmpzAlg.make) $@: No rule for such target" >&2
else
.DEFAULT::
	$(error PEDANTIC: $@: No rule for such target)
endif

PppmpzAlgclean ::
#-- end of cleanup_header ---------------
#-- start of cleanup_library -------------
	$(cleanup_echo) library PppmpzAlg
	-$(cleanup_silent) cd $(bin); /bin/rm -f $(library_prefix)PppmpzAlg$(library_suffix).a $(library_prefix)PppmpzAlg$(library_suffix).s? PppmpzAlg.stamp PppmpzAlg.shstamp
#-- end of cleanup_library ---------------
