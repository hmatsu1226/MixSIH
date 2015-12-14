#!/usr/bin/ruby

fin_prof = open(ARGV[0], "r")
fin_hap = open(ARGV[1], "r")
fout_prof = open(ARGV[2], "w")
fout_hap = open(ARGV[3], "w")
thresh = ARGV[4].to_f

tmp_prof = Array.new()
tmp_hap = Array.new()

fin_prof.gets
fin_hap.gets
phased = 0
while line1 = fin_prof.gets
	line2 = fin_hap.gets

	if line1[0,1] == '*'
		if tmp_prof.size != 0
			header = "BLOCK: offset: #{tmp_prof[0].split("\t")[0]} len: #{tmp_prof.size} phased: #{phased}"
			fout_prof.puts(header)
			fout_hap.puts(header)
		end

		for i in tmp_prof
			fout_prof.puts(i)
		end
		for i in tmp_hap
			fout_hap.puts(i)
		end
		if tmp_prof.size != 0
			fout_prof.puts("********")
			fout_hap.puts("********")
		end
		phased = 0
		tmp_prof.clear
		tmp_hap.clear
		fin_prof.gets
		fin_hap.gets

	elsif line1.split("\t")[3].to_f >= thresh
		if line2.split("\t")[1] != '-'
			phased += 1
		end
		tmp_prof.push(line1)
		tmp_hap.push(line2)

	else
		if tmp_prof.size != 0
			header = "BLOCK: offset: #{tmp_prof[0].split("\t")[0]} len: #{tmp_prof.size} phased: #{phased}"
			fout_prof.puts(header)
			fout_hap.puts(header)
		end

		for i in tmp_prof
			fout_prof.puts(i)
		end
		for i in tmp_hap
			fout_hap.puts(i)
		end
		if tmp_prof.size != 0
			fout_prof.puts("********")
			fout_hap.puts("********")
		end
		phased = 0
		tmp_prof.clear
		tmp_hap.clear

		if line2.split("\t")[1] != '-'
			phased += 1
		end

		tmp_prof.push(line1)
		tmp_hap.push(line2)
	end

end
