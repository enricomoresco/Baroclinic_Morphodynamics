# Changes with respect to ROMS

The configuration uses the unmodified ROMS source code (revision `57aecf5` of
<https://github.com/myroms/roms>). The only ROMS files that differ are the analytical functions in
`Functionals/`, which ROMS picks up instead of its own because `build_roms.sh` points
`MY_ANALYTICAL_DIR` to this folder, and the build script itself.

| File | Change |
|---|---|
| `ana_grid.h` | bed of the uniform flow with slope `user(6)`, plus the local (`user(7)`, sech²) and cumulative (`user(8)`, tanh) corrections |
| `ana_initial.h` | initial free surface with the same slope and cumulative rise; initial salinity equal to the tanh profile |
| `ana_fsobc.h` | no tide: sea level zero at the sea boundary |
| `ana_m2obc.h` | river discharge `user(1)*user(2)` per unit width at both ends, no tidal flow |
| `ana_tobc.h` | boundary salinity from the tanh profile |
| `ana_tclima.h` | salinity climatology (nudging target) equal to the tanh profile |
| `ana_sediment.h` | sediment bed 30 m thick with porosity 0.4, so that it cannot be depleted |
| `build_roms.sh` | application `ESTUARY_TEST`, serial gfortran build, no debugging, no PIO, `MY_ANALYTICAL_DIR` set to `Functionals/`, error if `ROMS_ROOT_DIR` is not set |

The `user(...)` parameters are listed in `estuary_test.h` and set by `scripts/make_run.py`.

## Unified diff

```diff
--- ROMS/Functionals/ana_grid.h
+++ roms/Functionals/ana_grid.h
@@ -954,7 +954,10 @@
 #elif defined ESTUARY_TEST
       DO j=JstrT,JendT
         DO i=IstrT,IendT
-          h(i,j)=5.0_r8+(Xsize-xr(i,j))/Xsize*5.0_r8
+          h(i,j)=user(1)-user(6)*xr(i,j)-                                   &
+     &           user(7)/COSH((xr(i,j)-user(4))/user(5))**2-            &
+     &           user(8)*(TANH((xr(i,j)-user(4))/user(5))-              &
+     &                    TANH(-user(4)/user(5)))
         END DO
       END DO
 #elif defined LAB_CANYON
--- ROMS/Functionals/ana_initial.h
+++ roms/Functionals/ana_initial.h
@@ -354,7 +354,15 @@
 !  Initial conditions for free-surface (m).
 !-----------------------------------------------------------------------
 !
-#if defined CHANNEL && !defined ONLY_TS_IC
+#if defined ESTUARY_TEST
+      DO j=JstrT,JendT
+        DO i=IstrT,IendT
+          zeta(i,j,1)=user(6)*xr(i,j)+                                   &
+     &                user(8)*(TANH((xr(i,j)-user(4))/user(5))-         &
+     &                         TANH(-user(4)/user(5)))
+        END DO
+      END DO
+#elif defined CHANNEL && !defined ONLY_TS_IC
       y0=0.5_r8*el(ng)
 # ifdef SOLVE3D
       DO j=JstrT,JendT
@@ -658,13 +666,7 @@
           DO i=IstrT,IendT
             t(i,j,k,1,itemp)=10.0_r8
 #  ifdef SALINITY
-            IF (xr(i,j).le.30000.0_r8) then
-              t(i,j,k,1,isalt)=30.0_r8
-            ELSEIF (xr(i,j).le.80000.0_r8) then
-              t(i,j,k,1,isalt)=(80000.0_r8-xr(i,j))/50000.0_r8*30.0_r8
-            ELSE
-              t(i,j,k,1,isalt)=0.0_r8
-            END IF
+            t(i,j,k,1,isalt)=0.5_r8*user(3)*(1.0_r8-TANH((xr(i,j)-user(4))/user(5)))
 #  endif
           END DO
         END DO
--- ROMS/Functionals/ana_fsobc.h
+++ roms/Functionals/ana_fsobc.h
@@ -106,7 +106,7 @@
 #elif defined ESTUARY_TEST
       IF (LBC(iwest,isFsur,ng)%acquire.and.                             &
      &    DOMAIN(ng)%Western_Edge(tile)) THEN
-        cff=1.0_r8*SIN(2.0_r8*pi*time(ng)/(12.0_r8*3600.0_r8))
+        cff=0.0_r8
         DO j=JstrT,JendT
           BOUNDARY(ng)%zeta_west(j)=cff
         END DO
--- ROMS/Functionals/ana_m2obc.h
+++ roms/Functionals/ana_m2obc.h
@@ -130,8 +130,8 @@
      &    DOMAIN(ng)%Western_Edge(tile)) THEN
         cff1=0.40_r8                                          ! west end
         cff2=0.08_r8
-        riv_flow=cff2*300.0_r8*5.0_r8
-        tid_flow=cff1*300.0_r8*10.0_r8
+        riv_flow=user(2)*user(1)*300.0_r8
+        tid_flow=0.0_r8
         my_area=0.0_r8
         my_flux=0.0_r8
         DO j=JstrP,JendP
@@ -150,8 +150,7 @@
       IF (LBC(ieast,isUbar,ng)%acquire.and.                             &
      &    LBC(ieast,isVbar,ng)%acquire.and.                             &
      &    DOMAIN(ng)%Eastern_Edge(tile)) THEN
-        cff2=0.08_r8                                          ! east end
-        riv_flow=cff2*300.0_r8*5.0_r8
+        riv_flow=user(2)*user(1)*300.0_r8                     ! east end
         my_area=0.0_r8
         my_flux=0.0_r8
         DO j=JstrP,JendP
--- ROMS/Functionals/ana_tobc.h
+++ roms/Functionals/ana_tobc.h
@@ -102,7 +102,7 @@
           DO j=JstrT,JendT
             BOUNDARY(ng)%t_east(j,k,itemp)=T0(ng)
 # ifdef SALINITY
-            BOUNDARY(ng)%t_east(j,k,isalt)=0.0_r8
+            BOUNDARY(ng)%t_east(j,k,isalt)=0.5_r8*user(3)*(1.0_r8-TANH((100000.0_r8-user(4))/user(5)))
 # endif
 # ifdef SEDIMENT
             DO ised=1,NST
@@ -119,7 +119,7 @@
           DO j=JstrT,JendT
             BOUNDARY(ng)%t_west(j,k,itemp)=T0(ng)
 # ifdef SALINITY
-            BOUNDARY(ng)%t_west(j,k,isalt)=30.0_r8
+            BOUNDARY(ng)%t_west(j,k,isalt)=0.5_r8*user(3)*(1.0_r8-TANH((0.0_r8-user(4))/user(5)))
 # endif
 # ifdef SEDIMENT
             DO ised=1,NST
--- ROMS/Functionals/ana_tclima.h
+++ roms/Functionals/ana_tclima.h
@@ -97,6 +97,18 @@
             END DO
           END DO
         END DO
+#elif defined ESTUARY_TEST
+        DO k=1,N(ng)
+          DO j=JstrT,JendT
+            DO i=IstrT,IendT
+              CLIMA(ng)%tclm(i,j,k,itemp)=10.0_r8
+# ifdef SALINITY
+              CLIMA(ng)%tclm(i,j,k,isalt)=0.5_r8*user(3)*               &
+     &          (1.0_r8-TANH((GRID(ng)%xr(i,j)-user(4))/user(5)))
+# endif
+            END DO
+          END DO
+        END DO
 #else
         DO k=1,N(ng)
           DO j=JstrT,JendT
--- ROMS/Functionals/ana_sediment.h
+++ roms/Functionals/ana_sediment.h
@@ -284,8 +284,8 @@
 !
           DO k=1,Nbed
             bed(i,j,k,iaged)=time(ng)
-            bed(i,j,k,ithck)=0.001_r8
-            bed(i,j,k,iporo)=0.90_r8
+            bed(i,j,k,ithck)=30.0_r8
+            bed(i,j,k,iporo)=0.40_r8
             DO ised=1,NST
               bed_frac(i,j,k,ised)=1.0_r8/REAL(NST,r8)
             END DO
--- ROMS/Bin/build_roms.sh
+++ roms/build_roms.sh
@@ -158,7 +158,7 @@
 # determine the name of the ".h" header file with the application
 # CPP definitions.
 
-export   ROMS_APPLICATION=UPWELLING
+export   ROMS_APPLICATION=ESTUARY_TEST
 
 # Set a local environmental variable to define the path to the directories
 # where the ROMS source code is located (MY_ROOT_DIR), and this project's
@@ -170,7 +170,12 @@
 if [ -n "${ROMS_ROOT_DIR:+1}" ]; then
   export      MY_ROOT_DIR=${ROMS_ROOT_DIR}
 else
-  export      MY_ROOT_DIR=${HOME}/ocean/repository/git
+  echo ""
+  echo "Set ROMS_ROOT_DIR to the directory that contains the ROMS source code"
+  echo "as a subdirectory named 'roms' (git clone of github.com/myroms/roms,"
+  echo "revision 57aecf5), e.g.   export ROMS_ROOT_DIR=\${HOME}/src"
+  echo ""
+  exit 1
 fi
 
 export     MY_PROJECT_DIR=${PWD}
@@ -226,8 +231,8 @@
 # out. Any string value (including off) will evaluate to TRUE in
 # conditional if-statements.
 
- export           USE_MPI=on            # distributed-memory parallelism
- export        USE_MPIF90=on            # compile with mpif90 script
+#export           USE_MPI=on            # distributed-memory parallelism
+#export        USE_MPIF90=on            # compile with mpif90 script
 #export         which_MPI=intel         # compile with mpiifort library
 #export         which_MPI=mpich         # compile with MPICH library
 #export         which_MPI=mpich2        # compile with MPICH2 library
@@ -238,12 +243,12 @@
 #export        USE_OpenMP=on            # shared-memory parallelism
 
 #export              FORT=ifx
- export              FORT=ifort
-#export              FORT=gfortran
+#export              FORT=ifort
+ export              FORT=gfortran
 #export              FORT=pgi
 
 if [ $g_flags -eq 1 ]; then
- export         USE_DEBUG=on            # use Fortran debugging flags
+  : # no debug
 fi
 
  export         USE_LARGE=on            # activate 64-bit compilation
@@ -271,7 +276,7 @@
 #export   USE_PARALLEL_IO=on            # Parallel I/O with NetCDF-4/HDF5
 
 if [ $pio_lib -eq 1 ]; then
- export           USE_PIO=on            # Parallel I/O with PIO library
+  : # no pio
 fi
 
 # If any of the coupling component use the HDF5 Fortran API for primary
@@ -310,7 +315,7 @@
 
  export     MY_HEADER_DIR=${MY_PROJECT_DIR}
 
- export MY_ANALYTICAL_DIR=${MY_PROJECT_DIR}
+ export MY_ANALYTICAL_DIR=${MY_PROJECT_DIR}/Functionals
 
 # Put the binary to execute in the following directory.
 
```
