/*
<pre>
Maintainer: Ulrich Küttler
            kuettler@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/Members/kuettler
            089 - 289-15238
</pre>
*/

for (vi=0; vi<8; ++vi)
{
#ifdef FLUID3_IS_TERM1
  /* Konvektionsterm */
  eforce_(vi*4)     += timefacfac*(velint_(0)*conv_r_(0, 0, vi) + velint_(1)*conv_r_(0, 1, vi) + velint_(2)*conv_r_(0, 2, vi)) ;
  eforce_(vi*4 + 1) += timefacfac*(velint_(0)*conv_r_(1, 0, vi) + velint_(1)*conv_r_(1, 1, vi) + velint_(2)*conv_r_(1, 2, vi)) ;
  eforce_(vi*4 + 2) += timefacfac*(velint_(0)*conv_r_(2, 0, vi) + velint_(1)*conv_r_(2, 1, vi) + velint_(2)*conv_r_(2, 2, vi)) ;
#endif

#ifdef FLUID3_IS_TERM2
  /* Stabilisierung der Konvektion ( L_conv_u) */
  eforce_(vi*4)     += 2.0*ttimetauM*conv_c_(vi)*conv_old_(0) ;
  eforce_(vi*4 + 1) += 2.0*ttimetauM*conv_c_(vi)*conv_old_(1) ;
  eforce_(vi*4 + 2) += 2.0*ttimetauM*conv_c_(vi)*conv_old_(2) ;
#endif

#ifdef FLUID3_IS_TERM3
  /* Stabilisierung der Konvektion (-L_visc_u) */
  eforce_(vi*4)     += -2.0*nu_*ttimetauM*conv_c_(vi)*visc_old_(0) ;
  eforce_(vi*4 + 1) += -2.0*nu_*ttimetauM*conv_c_(vi)*visc_old_(1) ;
  eforce_(vi*4 + 2) += -2.0*nu_*ttimetauM*conv_c_(vi)*visc_old_(2) ;
#endif

#ifdef FLUID3_IS_TERM4
  /* Viskositätsterm */
#endif

#ifdef FLUID3_IS_TERM5
  /* Stabilisierung der Viskosität ( L_conv_u) */
  eforce_(vi*4)     += 2.0*nu_*ttimetauMp*(conv_old_(0)*viscs2_(0, 0, vi) + conv_old_(1)*viscs2_(0, 1, vi) + conv_old_(2)*viscs2_(0, 2, vi)) ;
  eforce_(vi*4 + 1) += 2.0*nu_*ttimetauMp*(conv_old_(0)*viscs2_(0, 1, vi) + conv_old_(1)*viscs2_(1, 1, vi) + conv_old_(2)*viscs2_(1, 2, vi)) ;
  eforce_(vi*4 + 2) += 2.0*nu_*ttimetauMp*(conv_old_(0)*viscs2_(0, 2, vi) + conv_old_(1)*viscs2_(1, 2, vi) + conv_old_(2)*viscs2_(2, 2, vi)) ;
#endif

#ifdef FLUID3_IS_TERM6
  /* Stabilisierung der Viskosität (-L_visc_u) */
#endif

#ifdef FLUID3_IS_TERM7
  /* Druckterm */
#endif

#ifdef FLUID3_IS_TERM8
  /* Divergenzfreiheit */
#endif

#ifdef FLUID3_IS_TERM9
  /* Kontinuitätsstabilisierung */
#endif

#ifdef FLUID3_IS_TERM10
  /* Massenterm */
#endif

#ifdef FLUID3_IS_TERM11
  /* Konvektionsstabilisierung */
  eforce_(vi*4)     += timetauM*conv_c_(vi)*velint_(0) ;
  eforce_(vi*4 + 1) += timetauM*conv_c_(vi)*velint_(1) ;
  eforce_(vi*4 + 2) += timetauM*conv_c_(vi)*velint_(2) ;
#endif

#ifdef FLUID3_IS_TERM12
  /* Viskositätsstabilisierung */
#endif

#ifdef FLUID3_IS_TERM13
  /* Quellterm der rechten Seite */
  eforce_(vi*4)     += fac*funct_(vi)*rhsint_(0) ;
  eforce_(vi*4 + 1) += fac*funct_(vi)*rhsint_(1) ;
  eforce_(vi*4 + 2) += fac*funct_(vi)*rhsint_(2) ;
#endif

#ifdef FLUID3_IS_TERM14
  /* Konvektionsstabilisierung */
#endif

#ifdef FLUID3_IS_TERM15
  /* Viskositätsstabilisierung */
  eforce_(vi*4)     += 2.0*nu_*timetauMp*(rhsint_(0)*viscs2_(0, 0, vi) + rhsint_(1)*viscs2_(0, 1, vi) + rhsint_(2)*viscs2_(0, 2, vi)) ;
  eforce_(vi*4 + 1) += 2.0*nu_*timetauMp*(rhsint_(0)*viscs2_(0, 1, vi) + rhsint_(1)*viscs2_(1, 1, vi) + rhsint_(2)*viscs2_(1, 2, vi)) ;
  eforce_(vi*4 + 2) += 2.0*nu_*timetauMp*(rhsint_(0)*viscs2_(0, 2, vi) + rhsint_(1)*viscs2_(1, 2, vi) + rhsint_(2)*viscs2_(2, 2, vi)) ;
#endif

}

for (vi=8; vi<iel; ++vi)
{
#ifdef FLUID3_IS_TERM1
  /* Konvektionsterm */
  eforce_(vi*3 + 8) += timefacfac*(velint_(0)*conv_r_(0, 0, vi) + velint_(1)*conv_r_(0, 1, vi) + velint_(2)*conv_r_(0, 2, vi)) ;
  eforce_(vi*3 + 9) += timefacfac*(velint_(0)*conv_r_(1, 0, vi) + velint_(1)*conv_r_(1, 1, vi) + velint_(2)*conv_r_(1, 2, vi)) ;
  eforce_(vi*3 + 10) += timefacfac*(velint_(0)*conv_r_(2, 0, vi) + velint_(1)*conv_r_(2, 1, vi) + velint_(2)*conv_r_(2, 2, vi)) ;
#endif

#ifdef FLUID3_IS_TERM2
  /* Stabilisierung der Konvektion ( L_conv_u) */
  eforce_(vi*3 + 8) += 2.0*ttimetauM*conv_c_(vi)*conv_old_(0) ;
  eforce_(vi*3 + 9) += 2.0*ttimetauM*conv_c_(vi)*conv_old_(1) ;
  eforce_(vi*3 + 10) += 2.0*ttimetauM*conv_c_(vi)*conv_old_(2) ;
#endif

#ifdef FLUID3_IS_TERM3
  /* Stabilisierung der Konvektion (-L_visc_u) */
  eforce_(vi*3 + 8) += -2.0*nu_*ttimetauM*conv_c_(vi)*visc_old_(0) ;
  eforce_(vi*3 + 9) += -2.0*nu_*ttimetauM*conv_c_(vi)*visc_old_(1) ;
  eforce_(vi*3 + 10) += -2.0*nu_*ttimetauM*conv_c_(vi)*visc_old_(2) ;
#endif

#ifdef FLUID3_IS_TERM4
  /* Viskositätsterm */
#endif

#ifdef FLUID3_IS_TERM5
  /* Stabilisierung der Viskosität ( L_conv_u) */
  eforce_(vi*3 + 8) += 2.0*nu_*ttimetauMp*(conv_old_(0)*viscs2_(0, 0, vi) + conv_old_(1)*viscs2_(0, 1, vi) + conv_old_(2)*viscs2_(0, 2, vi)) ;
  eforce_(vi*3 + 9) += 2.0*nu_*ttimetauMp*(conv_old_(0)*viscs2_(0, 1, vi) + conv_old_(1)*viscs2_(1, 1, vi) + conv_old_(2)*viscs2_(1, 2, vi)) ;
  eforce_(vi*3 + 10) += 2.0*nu_*ttimetauMp*(conv_old_(0)*viscs2_(0, 2, vi) + conv_old_(1)*viscs2_(1, 2, vi) + conv_old_(2)*viscs2_(2, 2, vi)) ;
#endif

#ifdef FLUID3_IS_TERM6
  /* Stabilisierung der Viskosität (-L_visc_u) */
#endif

#ifdef FLUID3_IS_TERM7
  /* Druckterm */
#endif

#ifdef FLUID3_IS_TERM8
  /* Divergenzfreiheit */
#endif

#ifdef FLUID3_IS_TERM9
  /* Kontinuitätsstabilisierung */
#endif

#ifdef FLUID3_IS_TERM10
  /* Massenterm */
#endif

#ifdef FLUID3_IS_TERM11
  /* Konvektionsstabilisierung */
  eforce_(vi*3 + 8) += timetauM*conv_c_(vi)*velint_(0) ;
  eforce_(vi*3 + 9) += timetauM*conv_c_(vi)*velint_(1) ;
  eforce_(vi*3 + 10) += timetauM*conv_c_(vi)*velint_(2) ;
#endif

#ifdef FLUID3_IS_TERM12
  /* Viskositätsstabilisierung */
#endif

#ifdef FLUID3_IS_TERM13
  /* Quellterm der rechten Seite */
  eforce_(vi*3 + 8) += fac*funct_(vi)*rhsint_(0) ;
  eforce_(vi*3 + 9) += fac*funct_(vi)*rhsint_(1) ;
  eforce_(vi*3 + 10) += fac*funct_(vi)*rhsint_(2) ;
#endif

#ifdef FLUID3_IS_TERM14
  /* Konvektionsstabilisierung */
#endif

#ifdef FLUID3_IS_TERM15
  /* Viskositätsstabilisierung */
  eforce_(vi*3 + 8) += 2.0*nu_*timetauMp*(rhsint_(0)*viscs2_(0, 0, vi) + rhsint_(1)*viscs2_(0, 1, vi) + rhsint_(2)*viscs2_(0, 2, vi)) ;
  eforce_(vi*3 + 9) += 2.0*nu_*timetauMp*(rhsint_(0)*viscs2_(0, 1, vi) + rhsint_(1)*viscs2_(1, 1, vi) + rhsint_(2)*viscs2_(1, 2, vi)) ;
  eforce_(vi*3 + 10) += 2.0*nu_*timetauMp*(rhsint_(0)*viscs2_(0, 2, vi) + rhsint_(1)*viscs2_(1, 2, vi) + rhsint_(2)*viscs2_(2, 2, vi)) ;
#endif

}
