module mlp_namelist_read
    use ani_module, only: ani_nml
    implicit none

    contains

    subroutine read_mlp_namelist(read_file)
        logical, intent(in) :: read_file
        integer :: ifind, unit_number
        character(len=128) :: filename
        logical :: mdin_mlp
        integer :: anitype, ani_natom, mlp_model, mlp_shake
        character(len=8192) :: animask
        integer :: mlp_confine
        real :: confine1, confine2, confinek

        ! Define namelist
        namelist /mlp/ animask, anitype, animask, ani_natom, mlp_model, mlp_shake, mlp_confine, confine1, confine2, confinek

        ! Set default values
        mlp_shake=0
        anitype = 0
        animask = ''
        ani_natom = 0
        ani_nml%anitype = anitype
        ani_nml%animask = animask
        ani_nml%ani_natom = ani_natom
        ani_nml%mlp_model = 0 ! Default value, ANI-2x model
        ani_nml%mlp_shake = mlp_shake
        ani_nml%mlp_confine = 0
        ani_nml%confine1 = 0.0
        ani_nml%confine2 = ani_nml%confine1 + 1.0
        ani_nml%confinek = 0.0

        if (read_file) then
            rewind 5

            call nmlsrc('mlp',5,ifind)
            if (ifind /= 0) mdin_mlp=.true.

            !Read qmmm namelist
            rewind 5
            if ( mdin_mlp) then
                read(5,nml=mlp)
            else
                write(6, '(1x,a,/)') 'Could not find ani namelist'
                call mexit(6,1)
            endif
        endif

        ani_nml%anitype = anitype
        ani_nml%animask = animask
        ani_nml%ani_natom = ani_natom
        ani_nml%mlp_model = mlp_model
        ani_nml%mlp_shake = mlp_shake
        ani_nml%mlp_confine = mlp_confine
        if (mlp_confine == 1) then
            ani_nml%confine1 = confine1
            ani_nml%confine2 = confine2
            ani_nml%confinek = confinek

            write(6, *) 'Set confine1 = ', ani_nml%confine1, ' confine2 = ', ani_nml%confine2, ' confinek = ', ani_nml%confinek
        endif   

    end subroutine read_mlp_namelist

end module mlp_namelist_read
