PRO construct_trees

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
loop_start = 499 ;0
loop_stop  = 465 ;N_ELEMENTS(nrsubf)-1
FOR j=loop_start, loop_stop,-1 DO BEGIN  ; 99, 199, 299, 399, 499, 999, 1499, N_ELEMENTS(nrsubf)-1

    i_sub=nrsubf(j)
    str=STRTRIM(string(j),2)

	print,'HALO: ',j, '	OF	', loop_stop
	
	; open file for writing
	; ---------------------

	OPENW, lun, '/home/moon/sfortune/spinevo/newtrees/halo_'+str+'.dat', /GET_LUN
	printf,lun,format='("Most Massive Progenitor+All other Progenitors, Masses are at snap of max stellar mass for next Progenitors and gas mass from particles, First Progenitor mass is 2 snaps before identificationi")'
	printf,lun,format='(12(A14,2x))','SNAP','I_SUB','I_TREE','FILE_NR','REDSHIFT','M_STARS','M_STAR_2','M_GAS','M_GAS_2','M_DM','M_DM_2','MMP';,'X','Y','Z','VX','VY','VZ','MEAN_*_AGE','MEAN_METAL','M_STAR_TR','M_GAS_TR','M_DM_TR'

	find_halos_in_trees,i_sub,'136',0,t_haloes,i_tree			; find halo in trees at snap 136 simply matching positions

	; loop over all first progenitors
	; -------------------------------

	count=0.

    WHILE(i_tree GE 0)  DO BEGIN

		snap_fp=t_haloes[i_tree].SnapNum
		file_nr_fp=t_haloes[i_tree].FileNr
		i_sub_fp=t_haloes[i_tree].SubhaloIndex

		str1=string(snap_fp,form='(i03)')
		str2=STRTRIM(STRING(file_nr_fp),2)
		sub_name='/HydroSims/Magneticum/Box4/uhr_test/groups_'+str1+'/sub_'+str1+'.'+str2	
		readnew,sub_name,h,'HEAD'
		z=h.redshift
		HUBBLE=h.hubbleparam

		print,snap_fp

		m_star_fp=t_haloes[i_tree].SubMassTab[4]
		m_gas_fp=t_haloes[i_tree].SubMassTab[0]
		m_dm_fp=t_haloes[i_tree].SubMassTab[1]

		trace_back_first_prog,t_haloes,i_tree,m_star_fp_peak,m_gas_fp_peak,m_dm_fp_peak

		; write quantities of first progenitor here
		; -----------------------------------------

		printf,lun,format='(12(A14,2x))',snap_fp,i_sub_fp,i_tree,file_nr_fp,z,FLOAT(m_star_fp),FLOAT(m_star_fp_peak),FLOAT(m_gas_fp),FLOAT(m_gas_fp_peak),FLOAT(m_dm_fp),FLOAT(m_dm_fp_peak),1.;,pos[0,i_sub_fp],pos[1,i_sub_fp],pos[2,i_sub_fp],vel[0,i_sub_fp],vel[1,i_sub_fp],vel[2,i_sub_fp],mean_stellar_age[i_sub_fp],mean_metal[i_sub_fp],FLOAT(m_star_fp_trace),FLOAT(m_gas_fp_trace),FLOAT(m_dm_fp_trace)

		; loop over all progenitors
		; -------------------------

		i_np=t_haloes[i_tree].NextProgenitor

		WHILE(i_np GE 0) DO BEGIN
			
			snap_np=t_haloes[i_np].SnapNum
			file_nr_np=t_haloes[i_np].FileNr
			i_sub_np=t_haloes[i_np].SubhaloIndex

                        IF(snap_np NE snap_fp) THEN snap_np=snap_fp

                        m_star_np=t_haloes[i_np].SubMassTab[4]
                        m_gas_np=t_haloes[i_np].SubMassTab[0]
                        m_dm_np=t_haloes[i_np].SubMassTab[1]

			IF(m_star_np GT 0.0140800) THEN BEGIN

				trace_back_prog,t_haloes,i_np,m_star_np_peak,m_gas_np_peak,m_dm_np_peak
				
				; write quantities of next progenitors here
				; -----------------------------------------

				printf,lun,format='(12(A14,2x))',snap_np,i_sub_np,i_np,file_nr_np,z,FLOAT(m_star_np),FLOAT(m_star_np_peak),FLOAT(m_gas_np),FLOAT(m_gas_np_peak),FLOAT(m_dm_np),FLOAT(m_dm_np_peak),0.;,pos[0,i_sub_np],pos[1,i_sub_np],pos[2,i_sub_np],vel[0,i_sub_np],vel[1,i_sub_np],vel[2,i_sub_np];,mean_stellar_age[i_sub_np],mean_metal[i_sub_np]

			ENDIF	

			i_np=t_haloes[i_np].NextProgenitor	; set i_np to the index of the next progenitor
			
		ENDWHILE                                                                     
		
		; check whether halo jumps unreasonably (probably not needed for the new trees)
		; silvio: implementing Klaus' suggestion to switch back to central halo after switching to small subhalo

		x_old=t_haloes[i_tree].pos

		i_tree=t_haloes[i_tree].FirstProgenitor

		x_new=t_haloes[i_tree].pos

		r_dif=sqrt((x_old[0]-x_new[0])^2+(x_old[1]-x_new[1])^2+(x_old[2]-x_new[2])^2)

		IF(r_dif GT 200.) THEN BEGIN
			IF(i_tree NE t_haloes[i_tree].FIRSTHALOINFOFGROUP) THEN BEGIN
         		print,'Switching to main of FoF !'
				i_tree = t_haloes[i_tree].FIRSTHALOINFOFGROUP
				x_new=t_haloes[i_tree].pos
				r_dif=sqrt((x_old[0]-x_new[0])^2+(x_old[1]-x_new[1])^2+(x_old[2]-x_new[2])^2)
				IF(r_dif GT 200.) THEN BEGIN
					print,'Tree is broken even after switch:	Offset > ',r_dif,':	TERMINATE!'
					BREAK
				ENDIF
			ENDIF ELSE BEGIN
				print,'Tree is broken. Not a switch case.	Offset > ',r_dif,':	TERMINATE!'
				BREAK
			ENDELSE
		ENDIF

	ENDWHILE                                                                                                 		

	; next part is just if one wants to trace beyond snap 36 so in the next file
	; --------------------------------------------------------------------------

;	IF(snap_fp EQ 36) THEN BEGIN

;		find_halos_in_trees,i_sub_fp,'036',t_haloes2,i_tree

;		i_tree=t_haloes[i_tree].FirstProgenitor

;		WHILE(i_tree GE 0)  DO BEGIN

;			snap_fp=t_haloes2[i_tree].SnapNum
;			file_nr_fp=t_haloes2[i_tree].FileNr
;			i_sub_fp=t_haloes2[i_tree].SubhaloIndex

;			str1=string(snap_fp,form='(i03)')
;			str2=STRTRIM(STRING(file_nr_fp),2)
;			sub_name='/HydroSims/Magneticum/Box4/uhr_test/groups_'+str1+'/sub_'+str1+'.'+str2	
;			readnew,sub_name,h,'HEAD'
;			z=h.redshift
;			HUBBLE=h.hubbleparam

;			m_star_fp=t_haloes2[i_tree].SubMassTab[4]
;			m_gas_fp=t_haloes2[i_tree].SubMassTab[0]
;			m_dm_fp=t_haloes2[i_tree].SubMassTab[1]

;			trace_back_first_prog,t_haloes2,i_tree,m_star_fp_peak,m_gas_fp_peak,m_dm_fp_peak

			; write quantities of first progenitor here
			; -----------------------------------------

;			IF(snap_fp NE 36) THEN printf,lun,format='(12(A14,2x))',snap_fp,i_sub_fp,i_tree,file_nr_fp,z,FLOAT(m_star_fp),FLOAT(m_star_fp_peak),FLOAT(m_gas_fp),FLOAT(m_gas_fp_peak),FLOAT(m_dm_fp),FLOAT(m_dm_fp_peak),1.;,pos[0,i_sub_fp],pos[1,i_sub_fp],pos[2,i_sub_fp],vel[0,i_sub_fp],vel[1,i_sub_fp],vel[2,i_sub_fp],mean_stellar_age[i_sub_fp],mean_metal[i_sub_fp],FLOAT(m_star_fp_trace),FLOAT(m_gas_fp_trace),FLOAT(m_dm_fp_trace)

			; loop over all progenitors
			; -------------------------

;			i_np=t_haloes2[i_tree].NextProgenitor

;			WHILE(i_np GE 0) DO BEGIN
				
;				snap_np=t_haloes2[i_np].SnapNum
;				file_nr_np=t_haloes2[i_np].FileNr
;				i_sub_np=t_haloes2[i_np].SubhaloIndex

;				IF(snap_np NE snap_fp) THEN snap_np=snap_fp

;				m_star_np=t_haloes2[i_np].SubMassTab[4]
;				m_gas_np=t_haloes2[i_np].SubMassTab[0]
;				m_dm_np=t_haloes2[i_np].SubMassTab[1]

;				IF(m_star_np GT 0.0) THEN BEGIN

;					trace_back_prog,t_haloes2,i_np,m_star_np_peak,m_gas_np_peak,m_dm_np_peak,z_max
					
					; write quantities of next progenitors here
					; -----------------------------------------

;					printf,lun,format='(12(A14,2x))',snap_np,i_sub_np,i_np,file_nr_np,z,FLOAT(m_star_np),FLOAT(m_star_np_peak),FLOAT(m_gas_np),FLOAT(m_gas_np_peak),FLOAT(m_dm_np),FLOAT(m_dm_np_peak),0.;,pos[0,i_sub_np],pos[1,i_sub_np],pos[2,i_sub_np],vel[0,i_sub_np],vel[1,i_sub_np],vel[2,i_sub_np];,mean_stellar_age[i_sub_np],mean_metal[i_sub_np]

;				ENDIF	

;				i_np=t_haloes2[i_np].NextProgenitor
				
;			ENDWHILE                                                                                      
		
;			x_old=t_haloes2[i_tree].pos

;			i_tree=t_haloes2[i_tree].FirstProgenitor

;			x_new=t_haloes2[i_tree].pos

;			r_dif=sqrt((x_old[0]-x_new[0])^2+(x_old[1]-x_new[1])^2+(x_old[2]-x_new[2])^2)

;			IF(ABS(r_dif) GT 200.) THEN BEGIN
;				print,'Tree is broken: TERMINATE!'
;				BREAK
;			ENDIF

;			print,'SNAP: ',snap_fp

;		ENDWHILE
		
;	ENDIF

	FREE_LUN, lun
	CLOSE, lun

ENDFOR                                                                                                                 
        
stop

END

; example how to match from the trees to the snap files to also use particle information which is only
; available at snapshots in snap_all.
; Here one needs to be careful about the file number since t_haloes[*].SubhaloIndex contains the index within
; file t_haloes[*].FileNr and not in the total list!

PRO test_single_tree

snap_all=[136,132,128,124,120,116,112,108,106,104,102,100,96,92,88,84,80,76,68,64,60,58,52,48,44,40,36,32,28,24,20,16,12]

snr=string(136,form='(i03)')
read_tree,snr,t,t_haloes=t_haloes,ifile=0,sub='/HydroSims/Magneticum/Box4/uhr_test'

i_fp = 727509	; just some random halo to trace

WHILE(i_fp GT 0) DO BEGIN

    del=snap_all-t_haloes[i_fp].SnapNum
	mini=MIN(ABS(del))

	snr=string(t_haloes[i_fp].SnapNum,form='(i03)')
	filenr=STRTRIM(STRING(t_haloes[i_fp].FileNr),2)
	i_sub=t_haloes[i_fp].SubhaloIndex

	sub_name='/HydroSims/Magneticum/Box4/uhr_test/groups_'+snr+'/sub_'+snr+'.'+filenr
	snap_name='/HydroSims/Magneticum/Box4/uhr_test/snapdir_'+snr+'/snap_'+snr
	readnew,sub_name,h,'HEAD'
	z=h.redshift
	hubble=h.hubbleparam

	m_s = 0.
	m_g = 0.
	
	IF(mini EQ 0) THEN BEGIN

		read_halo_sim,snap_name,sub_name,i_sub,100.,id_s,x_stars,v_stars,m_stars,age_stars,Zs,iM,x_gas,v_gas,m_gas,gas_temp,sfr_gas,Zs_gas,x_dm,v_dm,gas_def,/cut_sub
;		plot,x_stars[0,*],x_stars[2,*],psym=1,symsize=0.1

		rr = sqrt( x_stars[0,*]^2+x_stars[1,*]^2+x_stars[2,*]^2 )
		kk = WHERE(rr LT 5.)
		m_s = total(m_stars(kk))

;                rr = sqrt(x_gas[0,*]^2+x_gas[1,*]^2+x_gas[2,*]^2)
;                kk = WHERE(rr LT 5.)
		kk = WHERE(gas_temp LT 1e5)
                m_g = total(m_gas(kk))


	ENDIF

	readnew,sub_name,grnr,'GRNR'
	readnew,sub_name,fsub,'FSUB'

;	IF(fsub(grnr(i_sub)) EQ i_sub) THEN type='main' ELSE type='sub'

	print,t_haloes[i_fp].SnapNum,t_haloes[i_fp].SubMassTab[4]*1e10/hubble,m_s,t_haloes[i_fp].SubMassTab[0]*1e10/hubble,m_g

	i_fp=t_haloes[i_fp].FirstProgenitor

ENDWHILE

STOP

END
