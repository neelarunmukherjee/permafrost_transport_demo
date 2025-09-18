module Reaction_Sandbox_neelresp_class
! Just trying to use one -neelarun mukherjee
#include "petsc/finclude/petscsys.h"
  use petscsys

  use Reaction_Sandbox_Base_class
  use PFLOTRAN_Constants_module

  implicit none

  private

  type, public, &
    extends(reaction_sandbox_base_type) :: reaction_sandbox_neelresp_type
    ! Aqueous species
    PetscInt :: species_DOM1_id
    PetscInt :: species_O2aq_id
    PetscInt :: species_HCO3_id
    PetscInt :: species_H_id
    ! PetscInt :: species_Eaq_id
    ! PetscInt :: species_Faq_id
    ! Immobile species (e.g. biomass)
    ! PetscInt :: species_Xim_id
    ! PetscInt :: species_Yim_id
  contains
    procedure, public :: Setup => neelrespSetup
    procedure, public :: Evaluate => neelrespEvaluate
  end type reaction_sandbox_neelresp_type

  public :: neelrespCreate

contains

! ************************************************************************** !

function neelrespCreate()
  !
  ! Allocates neelresp reaction object.
  ! Author: Glenn Hammond
  ! Date: 12/03/15
  ! Now neel is playig with it..

  implicit none

  class(reaction_sandbox_neelresp_type), pointer :: neelrespCreate

  allocate(neelrespCreate)
  nullify(neelrespCreate%next)

end function neelrespCreate

! ************************************************************************** !

subroutine neelrespSetup(this,reaction,option)
  !
  ! Sets up the neelresp reaction with hardwired parameters
  !
  ! Author: Glenn Hammond
  ! Date: 12/03/15

  use Reaction_Aux_module, only : reaction_rt_type, GetPrimarySpeciesIDFromName
  use Reaction_Immobile_Aux_module, only : GetImmobileSpeciesIDFromName
  use Option_module

  implicit none

  class(reaction_sandbox_neelresp_type) :: this
  class(reaction_rt_type) :: reaction
  type(option_type) :: option

  character(len=MAXWORDLENGTH) :: word

  ! Aqueous species
  word = 'DOM1'
  this%species_DOM1_id = &
    GetPrimarySpeciesIDFromName(word,reaction,option)
  word = 'O2(aq)'
  this%species_O2aq_id = &
    GetPrimarySpeciesIDFromName(word,reaction,option)
  word = 'HCO3-'
  this%species_HCO3_id = &
    GetPrimarySpeciesIDFromName(word,reaction,option)
  word = 'H+'
  this%species_H_id = &
    GetPrimarySpeciesIDFromName(word,reaction,option)
  ! word = 'Eaq'
  ! this%species_Eaq_id = &
  !   GetPrimarySpeciesIDFromName(word,reaction,option)
  ! word = 'Faq'
  ! this%species_Faq_id = &
  !   GetPrimarySpeciesIDFromName(word,reaction,option)

  ! Immobile species
  ! word = 'Xim'
  ! this%species_Xim_id = &
  !   GetImmobileSpeciesIDFromName(word,reaction%immobile,option)
  ! word = 'Yim'
  ! this%species_Yim_id = &
  !   GetImmobileSpeciesIDFromName(word,reaction%immobile,option)


end subroutine neelrespSetup

! ************************************************************************** !

subroutine neelrespEvaluate(this,Residual,Jacobian,compute_derivative, &
                          rt_auxvar,global_auxvar,material_auxvar,reaction, &
                          option)
  !
  ! Evaluates reaction storing residual and/or Jacobian
  !
  ! Author: Glenn Hammond
  ! Date: 12/03/15
  !
  use Option_module
  use Reaction_Aux_module
  use Reactive_Transport_Aux_module
  use Global_Aux_module
  use Material_Aux_module

  implicit none

  class(reaction_sandbox_neelresp_type) :: this
  type(option_type) :: option
  class(reaction_rt_type) :: reaction
  PetscBool :: compute_derivative
  ! the following arrays must be declared after reaction
  PetscReal :: Residual(reaction%ncomp)
  PetscReal :: Jacobian(reaction%ncomp,reaction%ncomp)
  type(reactive_transport_auxvar_type) :: rt_auxvar
  type(global_auxvar_type) :: global_auxvar
  type(material_auxvar_type) :: material_auxvar

  PetscInt, parameter :: iphase = 1
  PetscReal :: volume                 ! m^3 bulk volume
  PetscReal :: porosity               ! m^3 pore / m^3 bulk volume
  PetscReal :: liquid_saturation      ! m^3 water / m^3 pore space
  PetscReal :: L_water                ! L water

  PetscReal :: DOM1, O2aq, HCO3, H  !, Eaq, Faq  ! mol/L water
  ! PetscReal :: Xim, Yim  ! mol/m^3 bulk volume
  PetscReal :: Rate
  PetscReal :: RateDOM1, RateO2aq, RateHCO3, RateH !, RateE, RateF, RateX, RateY  ! mol/sec
  PetscReal :: stoichDOM1, stoichO2aq, stoichHCO3, stoichH !, stoichE, stoichF
  ! PetscReal :: stoichX, stoichY
  PetscReal :: n                    ! [-] Hill constant
  PetscReal :: k, kr        ! units are problem specific
  PetscReal :: K_DOM1, K_O2aq ! [mol/L water]

  porosity = material_auxvar%porosity               ! [m^3 pore/m^3 bulk volume]
  liquid_saturation = global_auxvar%sat(iphase)     ! [m^3 water/m^3 pore]
  volume = material_auxvar%volume                   ! [m^3 bulk volume]
            ! multiplying by 1.d3 converts m^3 water -> L water
  L_water = porosity*liquid_saturation*volume*1.d3

  DOM1 = rt_auxvar%total(this%species_DOM1_id,iphase) ! [mol/L water]
  O2aq = rt_auxvar%total(this%species_O2aq_id,iphase) ! [mol/L water]
  HCO3 = rt_auxvar%total(this%species_HCO3_id,iphase) ! [mol/L water]
  H = rt_auxvar%total(this%species_H_id,iphase) ! [mol/L water]
  ! Eaq = rt_auxvar%total(this%species_Eaq_id,iphase) ! [mol/L water]
  ! Faq = rt_auxvar%total(this%species_Faq_id,iphase) ! [mol/L water]

  ! Xim = rt_auxvar%immobile(this%species_Xim_id)     ! [mol/m^3 bulk volume]
  ! Yim = rt_auxvar%immobile(this%species_Yim_id)     ! [mol/m^3 bulk volume]

  ! initialize all rates to zero
  Rate = 0.d0

  RateDOM1 = 0.d0
  RateO2aq = 0.d0
  RateHCO3 = 0.d0
  stoichH = 0.d0
  ! RateE = 0.d0
  ! RateF = 0.d0
  ! RateX = 0.d0
  ! RateY = 0.d0

  ! stoichiometries
  ! reactants have negative stoichiometry
  ! products have positive stoichiometry
  stoichDOM1 = 0.d0
  stoichO2aq = 0.d0
  stoichHCO3 = 0.d0
  stoichH = 0.d0
  ! stoichE = 0.d0
  ! stoichF = 0.d0
  ! stoichX = 0.d0
  ! stoichY = 0.d0

  ! kinetic rate constants
  k = 0.d0
  kr = 0.d0

  ! Monod half-saturation constants
  K_DOM1 = 0.d0
  K_O2aq = 0.d0

  !----------------------------------------------------------------------------
  ! zero-order (A -> C)
  ! This rate constant is calculated to deplete a concentration of 1.d-3
  ! mol/liter in 25 years. The reaction must be run in batch mode (i.e. no
  ! transport). Increasing the rate constant will result in simulation failure
  ! unless the run time is reduced.
  !uncomment: k = 1.26839d-12  ! [mol/L water-sec]
  !uncomment: stoichA = -1.d0
  !uncomment: stoichC = 1.d0
  !uncomment: Rate = k * L_water  ! mol/sec
  !uncomment: RateA = stoichA * Rate
  !uncomment: RateC = stoichC * Rate

  !----------------------------------------------------------------------------
  ! first-order (A -> C)
  !uncomment: k = 1.d-9  ! [1/sec]
  !uncomment: stoichA = -1.d0
  !uncomment: stoichC = 1.d0
  !uncomment: Rate = k * Aaq * L_water  ! [mol/sec]
  !uncomment: RateA = stoichA * Rate
  !uncomment: RateC = stoichC * Rate

  !----------------------------------------------------------------------------
  ! second-order (A + B -> C)
  !uncomment: k = 1.d-6  ! [L water/mol-sec]
  !uncomment: stoichA = -1.d0
  !uncomment: stoichB = -1.d0
  !uncomment: stoichC = 1.d0
  !uncomment: Rate = k * Aaq * Baq * L_water  ! [mol/sec]
  !uncomment: RateA = stoichA * Rate
  !uncomment: RateB = stoichB * Rate
  !uncomment: RateC = stoichC * Rate

  !----------------------------------------------------------------------------
  ! Monod (A -> C)
  !uncomment: k = 1.d-12  ! [mol/L water-sec]
  !uncomment: K_Aaq = 1.d-4  ! [mol/L water]
  !uncomment: stoichA = -1.d0
  !uncomment: stoichC = 1.d0
  !uncomment: Rate = k * Aaq / (K_Aaq + Aaq) * L_water  ! [mol/sec]
  !uncomment: RateA = stoichA * Rate
  !uncomment: RateC = stoichC * Rate

  !----------------------------------------------------------------------------
  ! multiplicative Monod w/biomass  (A + B -> C + D)
  k = 1.d-6  ! [mol /mol biomass-sec]
  K_DOM1 = 1.d-4  ! [mol/L water]
  K_O2aq = 5.d-5  ! [mol/L water]
  stoichDOM1 = -1.d0
  stoichO2aq = -1.d0
  stoichHCO3 = 1.d0
  stoichH = 1.d0
  n = 1.2d0
  Rate = k  * DOM1 **n / (K_DOM1**n + DOM1**n) * &
                   O2aq / (K_O2aq + O2aq) * volume ! [mol/sec]
  RateDOM1 = stoichDOM1 * Rate
  RateO2aq = stoichO2aq * Rate
  RateHCO3 = stoichHCO3 * Rate
  RateH = stoichH * Rate

  !----------------------------------------------------------------------------
  ! decay and ingrowth (A -> B -> C)
  ! first-order rate constants
  !uncomment: k = 1.d-9  ! [1/sec]
  !uncomment: k1 = 2.d-9   ! [1/sec]
  !uncomment: k2 = 1.d-9  ! [1/sec]
  ! stoichiometries are moles of products generated from a mole of reactant
  !uncomment: stoichB = 3.d0
  !uncomment: stoichC = 1.d0
  !uncomment: stoichD = 2.d0
  !uncomment: Rate  = k *  Aaq * L_water  ! [mol/sec]
  !uncomment: Rate1 = k1 * Baq * L_water  ! [mol/sec]
  !uncomment: Rate2 = k2 * Caq * L_water  ! [mol/sec]
  !uncomment: RateA =                   - stoichB * Rate
  !uncomment: RateB =   stoichB * Rate  - stoichC * Rate1
  !uncomment: RateC =   stoichC * Rate1 - stoichD * Rate2
  !uncomment: RateD =   stoichD * Rate2

  !----------------------------------------------------------------------------
  ! first-order forward - reverse (A <-> C)
  !uncomment: k = 5.d-9   ! [1/sec]
  !uncomment: kr = 2.5d-9  ! [1/sec]
  !uncomment: stoichA = -1.d0
  !uncomment: stoichC = 1.d0
  !uncomment: Rate = (k * Aaq - kr * Caq) * L_water  ! [mol/sec]
  !uncomment: RateA = stoichA * Rate
  !uncomment: RateC = stoichC * Rate

  !----------------------------------------------------------------------------
  ! mass transfer between aqueous and immobile phases
  ! k [1/sec]
  ! kr [1/sec]
  ! Baq [mol/L water]
  ! Yim [mol/m^3 bulk volume]
  ! volume [m^3 bulk volume]
  ! L_water [L water]
  ! Rate [mole/sec]
  !uncomment: k = 5.d-10
  !uncomment: kr = 5.d-9
  !uncomment: stoichB = -5.d0
  !uncomment: stoichY = 1.d0
  !uncomment: Rate = k * Baq * L_water - kr * Yim * volume  ! [mol/sec]
  !uncomment: RateB = stoichB * Rate
  !uncomment: RateY = stoichY * Rate

  ! NOTES
  ! 1. Always subtract contribution from residual
  ! 2. Units of residual are moles/second
  Residual(this%species_DOM1_id) = Residual(this%species_DOM1_id) - RateDOM1
  Residual(this%species_O2aq_id) = Residual(this%species_O2aq_id) - RateO2aq
  Residual(this%species_HCO3_id) = Residual(this%species_HCO3_id) - RateHCO3
  Residual(this%species_H_id) = Residual(this%species_H_id) - RateH
  ! Residual(this%species_Eaq_id) = Residual(this%species_Eaq_id) - RateE
  ! Residual(this%species_Faq_id) = Residual(this%species_Faq_id) - RateF
  ! Residual(this%species_Xim_id + reaction%offset_immobile) = &
  !   Residual(this%species_Xim_id + reaction%offset_immobile) - RateX
  ! Residual(this%species_Yim_id + reaction%offset_immobile) = &
  !   Residual(this%species_Yim_id + reaction%offset_immobile) - RateY

  ! if (compute_derivative) then

  !   ! option%io_buffer = 'Reaction Sandbox neelresp does not support analytical &
  !   !                    &derivatives.'
  !   ! call PrintErrMsg(option)
  ! endif

end subroutine neelrespEvaluate

end module Reaction_Sandbox_neelresp_class
