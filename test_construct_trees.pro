
PRO test_construct_trees

; read trees
; ----------

snr=string(136,form='(i03)')
read_tree,snr,t,t_haloes=t_haloes,ifile=0,sub='/HydroSims/Magneticum/Box4/uhr_test' ; contains tree from snap 136 -> 36

;snr=string(036,form='(i03)')
;read_tree,snr,t,t_haloes=t_haloes2,ifile=0,sub='/HydroSims/Magneticum/Box4/uhr_test' ; contains tree from snap 36 -> 0
; (maybe the new trees are only one fiel)

file='/home/moon/fschulze/MA/Data/thesis/magneticum_data/rhea_list_snap_136_20180503.dat'
READCOLRHEA,file,nrsub,nrsubf,r_half_3d,r_half_2d,eps,eps_edge_on,l_r,l_r_edge_on,b,f_gas,triax,alph_jj,alph_mm,alph_jm,bet_jm_stars,bet_jm_gas,jx,jy,jz,m_star,m_tot,x,y,z,vx,vy,vz,x_hal,y_hal,z_hal,vx_hal,vy_hal,vz_hal,type,skipline=1,format='i,i,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,a'

; loop over all followed halos
; ----------------------------
loop_start = 2063 ;0
loop_stop  = 0 ;N_ELEMENTS(nrsubf)-1
FOR j=loop_start, loop_stop,-1 DO BEGIN  ; 99, 199, 299, 399, 499, 999, 1499, N_ELEMENTS(nrsubf)-1, for backwards add ',-1' after loop_stop

    i_sub=nrsubf(j)
    str=STRTRIM(string(j),2)

	find_halos_in_trees,i_sub,'136',0,t_haloes,i_tree	
	print,'HALO:',j, ' OF', loop_stop, ' I_SUB =', i_sub
	; sfc: central-only filter
	IF(i_tree NE t_haloes[i_tree].FIRSTHALOINFOFGROUP) THEN BEGIN
    	print,'Not a Central: TERMINATE'
	ENDIF ELSE BEGIN
    	print,'Central: GO ON'
		;; open file for writing
		;; ---------------------

		;;OPENW, lun, '/home/moon/sfortune/spinevo/test/halo_'+str+'.dat', /GET_LUN
		;;printf,lun,format='("Most Massive Progenitor+All other Progenitors, Masses are at snap of max stellar mass for next Progenitors and gas mass from particles, First Progenitor mass is 2 snaps before identificationi")'
		;;printf,lun,format='(12(A14,2x))','SNAP','I_SUB','I_TREE','FILE_NR','REDSHIFT','M_STARS','M_STAR_2','M_GAS','M_GAS_2','M_DM','M_DM_2','MMP';,'X','Y','Z','VX','VY','VZ','MEAN_*_AGE','MEAN_METAL','M_STAR_TR','M_GAS_TR','M_DM_TR'

		find_halos_in_trees,i_sub,'136',0,t_haloes,i_tree			; find halo in trees at snap 136 simply matching positions



		; loop over all first progenitors
		; ------------------------------

		count=0.
		fp_switch = 0
		fp_bordercross = 0

    	WHILE(i_tree GE 0)  DO BEGIN


			snap_fp=t_haloes[i_tree].SnapNum
			snap_old = snap_fp ; save to detect loop end later
			file_nr_fp=t_haloes[i_tree].FileNr
			i_sub_fp=t_haloes[i_tree].SubhaloIndex

			str1=string(snap_fp,form='(i03)')
			str2=STRTRIM(STRING(file_nr_fp),2)
			sub_name='/HydroSims/Magneticum/Box4/uhr_test/groups_'+str1+'/sub_'+str1+'.'+str2	
			readnew,sub_name,h,'HEAD'
			z=h.redshift
			HUBBLE=h.hubbleparam
			
			trace_back_first_prog,t_haloes,i_tree,m_star_fp_peak,m_gas_fp_peak,m_dm_fp_peak

			; Stay with central halos
			v_old = t_haloes[i_tree].vel * (1+z)
			x_old = t_haloes[i_tree].pos
			print,snap_fp,' SW =', fp_switch,' BW =', fp_bordercross

			i_tree=t_haloes[i_tree].FirstProgenitor
			snap_new = t_haloes[i_tree].SnapNum



			; fix for central halos and check for end with snap
			IF(i_tree NE t_haloes[i_tree].FIRSTHALOINFOFGROUP) THEN BEGIN
    	    	print,'Switching to main of FoF !'
				fp_switch = 1
				i_tree = t_haloes[i_tree].FIRSTHALOINFOFGROUP
				IF(snap_new GT snap_old) THEN BEGIN 
					print,'Snap Jump to ', snap_new, '	from ', snap_old
					BREAK
				ENDIF
			ENDIF ELSE BEGIN
				fp_switch = 0
			ENDELSE

			;;next snap
			;next_str1	= string(snap_new,form='(i03)')
			;next_fnr	= t_haloes[i_tree].FileNr
			;next_str2	= STRTRIM(STRING(next_fnr),2)
			;next_sub_name ='/HydroSims/Magneticum/Box4/uhr_test/groups_'+next_str1+'/sub_'+next_str1+'.'+next_str2	
			;readnew,next_sub_name,next_h,'HEAD'
			;next_z		= next_h.redshift
			;next_HUBBLE	= next_h.hubbleparam

			x_pred 	= x_old - 0.3*v_old			; 0.3 for maximum time-step
			x_new	= t_haloes[i_tree].pos
			r_pred	= sqrt((x_pred[0]-x_new[0])^2+(x_pred[1]-x_new[1])^2+(x_pred[2]-x_new[2])^2)
			r_dif	= sqrt((x_old[0]-x_new[0])^2+(x_old[1]-x_new[1])^2+(x_old[2]-x_new[2])^2)
			IF (r_dif-(0.3*NORM(v_old)) GT 0.) AND (r_dif LT 47500.) THEN BEGIN 
				print,'Radial codition: ', r_dif-(0.3*NORM(v_old)),' VEL =', NORM(v_old)
			ENDIF
			IF (r_pred GT 50) AND (r_dif LT 47500.) THEN BEGIN 
				print,'Prediction offset: ', r_pred,' VEL =', NORM(v_old)
			ENDIF
			IF (r_dif GT 47500.) THEN BEGIN 
				print,'Borderwalker: r_dif = ', r_dif
				fp_bordercross = 1
			ENDIF ELSE BEGIN
				fp_bordercross = 0
			ENDELSE
			IF (r_dif GT 200.) AND (r_dif LT 47500.) THEN BEGIN 
				print,'Distance Jump too big: TERMINATE! r_dif = ', r_dif
			ENDIF
		ENDWHILE
	ENDELSE

ENDFOR                                                                                                                 
        
stop

END
