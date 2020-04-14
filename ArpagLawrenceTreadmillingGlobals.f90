MODULE ArpagLawrenceTreadmillingGlobals
	
	IMPLICIT NONE
    INTEGER, PARAMETER  :: b8 = selected_real_kind(14) 

TYPE microtubule
    REAL(8) ::plusend
    REAL(8) ::minusend
    REAL(8) ::mtlength
    REAL(8) ::catstartp
    REAL(8) ::catstartm
    LOGICAL:: growingp
    LOGICAL:: shrinkingp
    LOGICAL:: growingm
    LOGICAL:: shrinkingm
    INTEGER:: dimercount
    INTEGER, dimension (:,:), ALLOCATABLE::nucleotide
    REAL(8) ::growthstartlp
    REAL(8) ::growthstarttp
    REAL(8) ::shrinkagestartlp
    REAL(8) ::shrinkagestarttp
    REAL(8) ::growthstartlm
    REAL(8) ::growthstarttm
    REAL(8) ::shrinkagestartlm
    REAL(8) ::shrinkagestarttm
    REAL(8) ::lengthinitial
    REAL(8) ::plusinitial
    REAL(8) ::minusinitial
    REAL(8) ::theta
    REAL(8) ::minusx
    REAL(8) ::minusy
    REAL(8) ::plusx
    REAL(8) ::plusy
    REAL(8) ::vgp
    REAL(8) ::vsp
    REAL(8) ::fcatp
    REAL(8) ::fresp
    REAL(8) ::vgm
    REAL(8) ::vsm
    REAL(8) ::fcatm
    REAL(8) ::fresm
    REAL(8) ::plusgrowthtime=0.0d0
    REAL(8) ::plusgrowthlength=0.0d0
    REAL(8) ::plusshrinkagetime=0.0d0
    REAL(8) ::plusshrinkagelength=0.0d0
    REAL(8) ::minusgrowthtime=0.0d0
    REAL(8) ::minusgrowthlength=0.0d0
    REAL(8) ::minusshrinkagetime=0.0d0
    REAL(8) ::minusshrinkagelength=0.0d0


END TYPE microtubule


END MODULE ArpagLawrenceTreadmillingGlobals
