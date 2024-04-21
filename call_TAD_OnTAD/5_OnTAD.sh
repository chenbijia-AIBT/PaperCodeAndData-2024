###############################5kb####################################################################
for r in 5000
do
	for i in GM12878 HMEC HUVEC IMR90 K562 NHEK
	do
		mkdir -p /mnt/f/project1/12_revise/TAD_calling/methods/0.OnTAD/1.OnTAD_output/${i}/${r}
		chr_int=/mnt/h/bxm/Ref/chrom_hg19.sizes
		matrix_dr=/mnt/f/project1/12_revise/TAD_calling/0.contact_matrices/Input/${i}/ICE/dense/${r}/
		out_dr=/mnt/f/project1/12_revise/TAD_calling/methods/0.OnTAD/1.OnTAD_output/${i}/${r}/
		cd ${out_dr}
		#遍历文件夹matrix_dr中的所有文件
		for ii in $( ls ${matrix_dr})
		do
			#定义chr
			chr=${ii%%_*}
			#获取每条染色体对应的chromosome_size
			chr_size=$(cat ${chr_int} | awk -v OFS="," '{print $1,$2}' | grep ${chr}, | awk -v FS=',' '{print $2}')
			echo ${ii} is ${chr},size is ${chr_size}
			# ##的意思是截掉字符串chr
			/mnt/f/project1/1_OnTAD/1-OnTAD/OnTAD-master/OnTAD ${matrix_dr}${ii} -penalty 0.1 -maxsz 400 -minsz 6 -o OnTAD_ICE_pen0.1_max_400_min_6_${chr} 
		done 
	done
done
###############################10kb####################################################################
for r in 10000
do
	for i in GM12878 HMEC HUVEC IMR90 K562 NHEK
	do
		mkdir -p /mnt/f/project1/12_revise/TAD_calling/methods/0.OnTAD/1.OnTAD_output/${i}/${r}
		chr_int=/mnt/h/bxm/Ref/chrom_hg19.sizes
		matrix_dr=/mnt/f/project1/12_revise/TAD_calling/0.contact_matrices/Input/${i}/ICE/dense/${r}/
		out_dr=/mnt/f/project1/12_revise/TAD_calling/methods/0.OnTAD/1.OnTAD_output/${i}/${r}/
		cd ${out_dr}
		#遍历文件夹matrix_dr中的所有文件
		for ii in $( ls ${matrix_dr})
		do
			#定义chr
			chr=${ii%%_*}
			#获取每条染色体对应的chromosome_size
			chr_size=$(cat ${chr_int} | awk -v OFS="," '{print $1,$2}' | grep ${chr}, | awk -v FS=',' '{print $2}')
			echo ${ii} is ${chr},size is ${chr_size}
			# ##的意思是截掉字符串chr
			/mnt/f/project1/1_OnTAD/1-OnTAD/OnTAD-master/OnTAD ${matrix_dr}${ii} -penalty 0.1 -maxsz 200 -minsz 3 -o OnTAD_ICE_pen0.1_max_200_min_3_${chr} 
		done 
	done
done
###############################25kb####################################################################
for r in 25000
do
	for i in GM12878 HMEC HUVEC IMR90 K562 NHEK
	do
		mkdir -p /mnt/f/project1/12_revise/TAD_calling/methods/0.OnTAD/1.OnTAD_output/${i}/${r}
		chr_int=/mnt/h/bxm/Ref/chrom_hg19.sizes
		matrix_dr=/mnt/f/project1/12_revise/TAD_calling/0.contact_matrices/Input/${i}/ICE/dense/${r}/
		out_dr=/mnt/f/project1/12_revise/TAD_calling/methods/0.OnTAD/1.OnTAD_output/${i}/${r}/
		cd ${out_dr}
		#遍历文件夹matrix_dr中的所有文件
		for ii in $( ls ${matrix_dr})
		do
			#定义chr
			chr=${ii%%_*}
			#获取每条染色体对应的chromosome_size
			chr_size=$(cat ${chr_int} | awk -v OFS="," '{print $1,$2}' | grep ${chr}, | awk -v FS=',' '{print $2}')
			echo ${ii} is ${chr},size is ${chr_size}
			# ##的意思是截掉字符串chr
			/mnt/f/project1/1_OnTAD/1-OnTAD/OnTAD-master/OnTAD ${matrix_dr}${ii} -penalty 0.1 -maxsz 80 -minsz 1.2 -o OnTAD_ICE_pen0.1_max_80_min_1.2_${chr} 
		done 
	done
done
###############################50kb####################################################################
for r in 50000
do
	for i in GM12878 HMEC HUVEC IMR90 K562 NHEK
	do
		mkdir -p /mnt/f/project1/12_revise/TAD_calling/methods/0.OnTAD/1.OnTAD_output/${i}/${r}
		chr_int=/mnt/h/bxm/Ref/chrom_hg19.sizes
		matrix_dr=/mnt/f/project1/12_revise/TAD_calling/0.contact_matrices/Input/${i}/ICE/dense/${r}/
		out_dr=/mnt/f/project1/12_revise/TAD_calling/methods/0.OnTAD/1.OnTAD_output/${i}/${r}/
		cd ${out_dr}
		#遍历文件夹matrix_dr中的所有文件
		for ii in $( ls ${matrix_dr})
		do
			#定义chr
			chr=${ii%%_*}
			#获取每条染色体对应的chromosome_size
			chr_size=$(cat ${chr_int} | awk -v OFS="," '{print $1,$2}' | grep ${chr}, | awk -v FS=',' '{print $2}')
			echo ${ii} is ${chr},size is ${chr_size}
			# ##的意思是截掉字符串chr
			/mnt/f/project1/1_OnTAD/1-OnTAD/OnTAD-master/OnTAD ${matrix_dr}${ii} -penalty 0.1 -maxsz 40 -minsz 0.6 -o OnTAD_ICE_pen0.1_max_40_min_0.6_${chr} 
		done 
	done
done


