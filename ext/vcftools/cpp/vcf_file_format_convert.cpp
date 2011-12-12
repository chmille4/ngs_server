/*
 * vcf_file_output.cpp
 *
 *  Created on: Aug 28, 2009
 *      Author: Adam Auton
 *      ($Revision: 249 $)
 */
#include "vcf_file.h"
/*
void vcf_file::output_as_plink(const string &output_file_prefix)
{
	// Output as PLINK formatted PED/MAP files.
	if (has_genotypes == false)
		error("Require Genotypes in VCF file in order to output as PLINK.");

	printLOG("Writing PLINK PED file ... ");
	string ped_file = output_file_prefix + ".ped";
	string map_file = output_file_prefix + ".map";

	ofstream PED(ped_file.c_str());
	if (!PED.is_open()) error("Could not open output file: " + ped_file, 12);

	vector<string> alleles;
	char phase;
	pair<int, int> genotype;
	string vcf_line;
	vcf_entry e(N_indv);
	for (unsigned int ui=0; ui<N_indv; ui++)
	{
		if (include_indv[ui] == false)
			continue;
		PED << indv[ui] << "\t" << indv[ui] << "\t" << 0 << "\t" << 0 << "\t" << 0 << "\t" << 0;

		for (unsigned int s=0; s<N_entries; s++)
		{
			if (include_entry[s] == false)
				continue;

			get_vcf_entry(s, vcf_line);
			e.reset(vcf_line);
			e.parse_basic_entry(true);

			if (e.get_N_alleles() <= 2)	// Only output sites with one alternative allele
			{
				e.get_alleles_vector(alleles);
				genotype = make_pair(-1,-1);
				phase = '/';
				if (include_genotype[s][ui] == true)
				{
					e.parse_genotype_entry(ui, true);
					e.get_indv_GENOTYPE_ids(ui, genotype);
					phase = e.get_indv_PHASE(ui);
				}

				if (genotype.first == -1)
					PED << "\t0";
				else
					PED << "\t" << alleles[genotype.first];

				if (genotype.second == -1)
				{
					if (phase == '/')
						PED << "\t0";
					else if (genotype.first != -1)
						PED << "\t" << alleles[genotype.first];	// Male X-chr, Y-chr etc
					else
						PED << "\t0";
				}
				else
					PED << "\t" << alleles[genotype.second];
			}
		}
		PED << endl;
	}

	PED.close();

	printLOG("Writing PLINK MAP file ... ");
	ofstream MAP(map_file.c_str());
	if (!MAP.is_open()) error("Could not open output file: " + map_file, 12);
	for (unsigned int s=0; s<N_entries; s++)
	{
		if (include_entry[s] == false)
			continue;

		get_vcf_entry(s, vcf_line);
		e.reset(vcf_line);
		e.parse_basic_entry(true);
		if (e.get_N_alleles() <= 2)	// Only output sites with one alternative allele
		{
			if (e.get_ID() == ".")
				MAP << e.get_CHROM() << "\t" << e.get_POS() << "\t0\t" << e.get_POS() << endl;
			else
				MAP << e.get_CHROM() << "\t" << e.get_ID() << "\t0\t" << e.get_POS() << endl;
		}
	}

	MAP.close();
	printLOG("Done.\n");
}
*/

void vcf_file::output_as_plink(const string &output_file_prefix)
{
	// Output as PLINK formatted PED/MAP files.
	if (has_genotypes == false)
		error("Require Genotypes in VCF file in order to output as PLINK.");

	printLOG("Writing PLINK PED file ... \n");
	string ped_file = output_file_prefix + ".ped";
	string map_file = output_file_prefix + ".map";

	vector<ofstream *> tmp_files(N_indv);
	vector<string> tmp_filenames(N_indv);
	for (unsigned int ui=0; ui<N_indv; ui++)
	{
		if (include_indv[ui] == false)
			continue;
		string filename(tmpnam(NULL));
		ofstream *tmp_file = new ofstream(filename.c_str());
		if (!tmp_file->good())
			error("\n\nCould not open temporary file.\n\n"
					"Most likely this is because the system is not allowing me to open enough temporary files.\n"
					"Try using ulimit -n <int> to increase the number of allowed open files.\n"
					"Alternatively, try the --plink-tped command.", 12);
		(*tmp_file) << indv[ui] << "\t" << indv[ui] << "\t" << 0 << "\t" << 0 << "\t" << 0 << "\t" << 0;
		tmp_files[ui] = tmp_file;
		tmp_filenames[ui] = filename;
	}

	vector<string> alleles;
	char phase;
	pair<int, int> genotype;
	string vcf_line;
	vcf_entry e(N_indv);
	ofstream *tmp_file;
	for (unsigned int s=0; s<N_entries; s++)
	{
		if (include_entry[s] == false)
			continue;

		get_vcf_entry(s, vcf_line);
		e.reset(vcf_line);
		e.parse_basic_entry(true);

		if (e.get_N_alleles() > 2)
		{
			one_off_warning("\tPLINK: Only outputting biallelic loci.");
			continue;
		}

		e.get_alleles_vector(alleles);

		for (unsigned int ui=0; ui<N_indv; ui++)
		{
			if (include_indv[ui] == false)
				continue;

			tmp_file = tmp_files[ui];

			genotype = make_pair(-1,-1);
			phase = '/';
			if (include_genotype[s][ui] == true)
			{
				e.parse_genotype_entry(ui, true);
				e.get_indv_GENOTYPE_ids(ui, genotype);
				phase = e.get_indv_PHASE(ui);
			}

			if (genotype.first == -1)
				(*tmp_file) << "\t0";
			else
				(*tmp_file) << "\t" << alleles[genotype.first];

			if (genotype.second == -1)
			{
				if (phase == '/')
					(*tmp_file) << "\t0";
				else if (genotype.first != -1)
					(*tmp_file) << "\t" << alleles[genotype.first];	// Male X-chr, Y-chr etc
				else
					(*tmp_file) << "\t0";
			}
			else
				(*tmp_file) << "\t" << alleles[genotype.second];
		}
	}

	ofstream PED(ped_file.c_str());
	if (!PED.is_open()) error("Could not open output file: " + ped_file, 12);
	string tmp_line;
	for (unsigned int ui=0; ui<N_indv; ui++)
	{
		if (include_indv[ui] == false)
			continue;

		ofstream *tmp_file = tmp_files[ui];
		(*tmp_file) << endl;
		tmp_file->close();

		ifstream read_file(tmp_filenames[ui].c_str());
		if (!read_file.good())
			error("\n\nCould not open temporary file.\n\n"
					"Most likely this is because the system is not allowing me to open enough temporary files.\n"
					"Try using ulimit -n <int> to increase the number of allowed open files.\n"
					"Alternatively, try the --plink-tped command.", 12);
		getline(read_file, tmp_line);
		PED << tmp_line << endl;
		read_file.close();
		remove(tmp_filenames[ui].c_str());
	}
	PED.close();

	printLOG("Writing PLINK MAP file ... ");
	ofstream MAP(map_file.c_str());
	if (!MAP.is_open()) error("Could not open output file: " + map_file, 12);
	int POS; string ID;
	for (unsigned int s=0; s<N_entries; s++)
	{
		if (include_entry[s] == false)
			continue;

		get_vcf_entry(s, vcf_line);
		e.reset(vcf_line);
		e.parse_basic_entry(true);
		if (e.get_N_alleles() > 2)
			continue;
		POS = e.get_POS();
		ID = e.get_ID();
		if (ID == ".")
			MAP << e.get_CHROM() << "\t" << POS << "\t0\t" << POS << endl;
		else
			MAP << e.get_CHROM() << "\t" << ID << "\t0\t" << POS << endl;
	}

	MAP.close();
	printLOG("Done.\n");
}

// Output as Plink Transposed file
void vcf_file::output_as_plink_tped(const string &output_file_prefix)
{
	// Output as PLINK formatted PED/MAP files.
	if (has_genotypes == false)
		error("Require Genotypes in VCF file in order to output as PLINK TPED.");

	printLOG("Writing PLINK TPED file ... ");
	string tped_file = output_file_prefix + ".tped";
	string tfam_file = output_file_prefix + ".tfam";

	ofstream TPED(tped_file.c_str());
	if (!TPED.is_open()) error("Could not open output file: " + tped_file, 12);

	vector<string> alleles;
	char phase;
	pair<int, int> genotype;
	string vcf_line;
	vcf_entry e(N_indv);
	for (unsigned int s=0; s<N_entries; s++)
	{
		if (include_entry[s] == false)
			continue;

		get_vcf_entry(s, vcf_line);
		e.reset(vcf_line);
		e.parse_basic_entry(true);

		if (e.get_N_alleles() > 2)	// Only output sites with at most one alternative allele
		{
			one_off_warning("\tPLINK-TPED: Only outputting biallelic loci.");
			continue;
		}

		if (e.get_ID() == ".")
			TPED << e.get_CHROM() << "\t" << e.get_POS() << "\t0\t" << e.get_POS();
		else
			TPED << e.get_CHROM() << "\t" << e.get_ID() << "\t0\t" << e.get_POS();

		e.get_alleles_vector(alleles);

		for (unsigned int ui=0; ui<N_indv; ui++)
		{
			if (include_indv[ui] == false)
				continue;

			genotype = make_pair(-1,-1);
			phase = '/';
			if (include_genotype[s][ui] == true)
			{
				e.parse_genotype_entry(ui, true);
				e.get_indv_GENOTYPE_ids(ui, genotype);
				phase = e.get_indv_PHASE(ui);
			}

			if (genotype.first == -1)
				TPED << "\t0";
			else
				TPED << "\t" << alleles[genotype.first];

			if (genotype.second == -1)
			{
				if (phase == '/')
					TPED << "\t0";
				else if (genotype.first != -1)
					TPED << "\t" << alleles[genotype.first];	// Male X-chr, Y-chr etc
				else
					TPED << "\t0";
			}
			else
				TPED << "\t" << alleles[genotype.second];
		}
		TPED << endl;
	}

	TPED.close();

	printLOG("Writing PLINK TFAM file ... ");
	ofstream TFAM(tfam_file.c_str());
	if (!TFAM.is_open()) error("Could not open output file: " + tfam_file, 12);
	for (unsigned int ui=0; ui<N_indv; ui++)
	{
		if (include_indv[ui] == false)
			continue;
		TFAM << indv[ui] << "\t" << indv[ui] << "\t" << 0 << "\t" << 0 << "\t" << 0 << "\t" << 0 << endl;
	}

	TFAM.close();
	printLOG("Done.\n");
}

/*
// Output as a simple 0/1/2 matrix
void vcf_file::output_as_012_matrix(const string &output_file_prefix)
{
	if (has_genotypes == false)
		error("Require Genotypes in VCF file in order to output as 0/1/2 matrix.");

	printLOG("Writing 012 matrix file ... ");
	string ped_file = output_file_prefix + ".012";
	string map_file = output_file_prefix + ".012.pos";
	string fam_file = output_file_prefix + ".012.indv";

	ofstream PED(ped_file.c_str());
	if (!PED.is_open()) error("Could not open output file: " + ped_file, 12);
	string allele1, allele2;

	ofstream FAM(fam_file.c_str());
	if (!FAM.is_open()) error("Could not open output file: " + fam_file, 12);

	pair<int, int> genotype;
	string vcf_line;
	vcf_entry e(N_indv);
	for (unsigned int ui=0; ui<N_indv; ui++)
	{
		if (include_indv[ui] == false)
			continue;
		FAM << indv[ui] << endl;
		PED << ui;
		//uk = 2*ui;
		for (unsigned int s=0; s<N_entries; s++)
		{
			if (include_entry[s] == false)
				continue;

			get_vcf_entry(s, vcf_line);
			e.reset(vcf_line);
			e.parse_basic_entry(true);

			if (e.get_N_alleles() <= 2)	// Only output sites with one alternative allele
			{
				genotype = make_pair(-1,-1);
				if (include_genotype[s][ui] == true)
				{
					e.parse_genotype_entry(ui, true);
					e.get_indv_GENOTYPE_ids(ui, genotype);
				}

				if ((genotype.first == -1) && (genotype.second == -1))
					PED << "\t-1";	// Missing data
				else if ((genotype.first == 0) && (genotype.second == 0))
					PED << "\t0";	// No copies of the alternative allele
				else
				{
					if ((genotype.first == 1) && (genotype.second == 1))
						PED << "\t2";	// Two copies of the alternative allele
					else
						PED << "\t1";	// Must be one copy of the alternative allele.
				}
			}
		}
		PED << endl;
	}

	FAM.close();
	PED.close();

	ofstream MAP(map_file.c_str());
	if (!MAP.is_open()) error("Could not open output file: " + map_file, 12);
	for (unsigned int s=0; s<N_entries; s++)
	{
		if (include_entry[s] == false)
			continue;

		get_vcf_entry(s, vcf_line);
		e.reset(vcf_line);
		e.parse_basic_entry(true);
		if (e.get_N_alleles() <= 2)	// Only output sites with one alternative allele
		{
			MAP << e.get_CHROM() << "\t" << e.get_POS() << endl;
		}
	}

	MAP.close();
	printLOG("Done.\n");
}
*/

void vcf_file::output_as_012_matrix(const string &output_file_prefix)
{
	// Output as PLINK formatted PED/MAP files.
	if (has_genotypes == false)
		error("Require Genotypes in VCF file in order to output as 0/1/2 matrix.");

	printLOG("Writing 012 matrix file ... ");
	string ped_file = output_file_prefix + ".012";
	string map_file = output_file_prefix + ".012.pos";
	string fam_file = output_file_prefix + ".012.indv";

	ofstream FAM(fam_file.c_str());
	if (!FAM.is_open()) error("Could not open output file: " + fam_file, 12);

	vector<ofstream *> tmp_files(N_indv);
	vector<string> tmp_filenames(N_indv);
	for (unsigned int ui=0; ui<N_indv; ui++)
	{
		if (include_indv[ui] == false)
			continue;
		FAM << indv[ui] << endl;
		string filename(tmpnam(NULL));
		ofstream *tmp_file = new ofstream(filename.c_str());
		if (!tmp_file->good())
			error("\n\nCould not open temporary file.\n\n"
			"Most likely this is because the system is not allowing me to open enough temporary files.\n"
			"Try using ulimit -n <int> to increase the number of allowed open files.\n", 12);
		(*tmp_file) << ui;
		tmp_files[ui] = tmp_file;
		tmp_filenames[ui] = filename;
	}

	FAM.close();

	vector<string> alleles;
	char phase;
	pair<int, int> genotype;
	string vcf_line;
	vcf_entry e(N_indv);
	ofstream *tmp_file;
	for (unsigned int s=0; s<N_entries; s++)
	{
		if (include_entry[s] == false)
			continue;

		get_vcf_entry(s, vcf_line);
		e.reset(vcf_line);
		e.parse_basic_entry(true);

		if (e.get_N_alleles() > 2)
		{
			one_off_warning("\t012: Only outputting biallelic loci.");
			continue;
		}

		e.get_alleles_vector(alleles);

		for (unsigned int ui=0; ui<N_indv; ui++)
		{
			if (include_indv[ui] == false)
				continue;

			tmp_file = tmp_files[ui];

			genotype = make_pair(-1,-1);
			phase = '/';
			if (include_genotype[s][ui] == true)
			{
				e.parse_genotype_entry(ui, true);
				e.get_indv_GENOTYPE_ids(ui, genotype);
				phase = e.get_indv_PHASE(ui);
			}

			if ((genotype.first == -1) && (genotype.second == -1))
				(*tmp_file) << "\t-1";	// Missing data
			else if ((genotype.first == 0) && (genotype.second == 0))
				(*tmp_file) << "\t0";	// No copies of the alternative allele
			else
			{
				if ((genotype.first == 1) && (genotype.second == 1))
					(*tmp_file) << "\t2";	// Two copies of the alternative allele
				else
					(*tmp_file) << "\t1";	// Must be one copy of the alternative allele.
			}
		}
	}

	ofstream PED(ped_file.c_str());
	if (!PED.is_open()) error("Could not open output file: " + ped_file, 12);
	string tmp_line;
	for (unsigned int ui=0; ui<N_indv; ui++)
	{
		if (include_indv[ui] == false)
			continue;

		ofstream *tmp_file = tmp_files[ui];
		(*tmp_file) << endl;
		tmp_file->close();

		ifstream read_file(tmp_filenames[ui].c_str());
		if (!read_file.good())
			error("\n\nCould not open temporary file.\n\n"
			"Most likely this is because the system is not allowing me to open enough temporary files.\n"
			"Try using ulimit -n <int> to increase the number of allowed open files.\n", 12);
		getline(read_file, tmp_line);
		PED << tmp_line << endl;
		read_file.close();
		remove(tmp_filenames[ui].c_str());
	}
	PED.close();

	printLOG("Writing 012 positions file ... ");
	ofstream MAP(map_file.c_str());
	if (!MAP.is_open()) error("Could not open output file: " + map_file, 12);
	for (unsigned int s=0; s<N_entries; s++)
	{
		if (include_entry[s] == false)
			continue;

		get_vcf_entry(s, vcf_line);
		e.reset(vcf_line);
		e.parse_basic_entry(true);
		if (e.get_N_alleles() <= 2)	// Only output sites with one alternative allele
		{
			MAP << e.get_CHROM() << "\t" << e.get_POS() << endl;
		}
	}

	MAP.close();
	printLOG("Done.\n");
}

// Output as IMPUTE format
void vcf_file::output_as_IMPUTE(const string &output_file_prefix)
{
	if (has_genotypes == false)
		error("Require Genotypes in VCF file in order to output IMPUTE format.");

	printLOG("Outputting in IMPUTE format (bi-allelic, completely phased SNPs only)\n");
	unsigned int s, ui;
	string legend_file = output_file_prefix + ".impute.legend";
	string haplotype_file = output_file_prefix + ".impute.hap";
	string indv_file = output_file_prefix + ".impute.hap.indv";
	ofstream legend(legend_file.c_str());
	if (!legend.is_open())
		error("Could not open IMPUTE Legend Output File: " + legend_file, 2);
	legend << "ID pos allele0 allele1" << endl;

	ofstream hap(haplotype_file.c_str());
	if (!hap.is_open())
		error("Could not open IMPUTE Haplotype Output File: " + haplotype_file, 2);

	ofstream indv_out(indv_file.c_str());
	if (!indv_out.is_open())
		error("Could not open IMPUTE Individual Output File: " + indv_file, 2);

	for (ui=0; ui<N_indv; ui++)
	{
		if (include_indv[ui] == false)
			continue;
		indv_out << indv[ui] << endl;
	}
	indv_out.close();

	pair<int, int> alleles;
	string vcf_line;
	vcf_entry e(N_indv);
	for (s=0; s<N_entries; s++)
	{
		if (include_entry[s] == false)
			continue;

		get_vcf_entry(s, vcf_line);
		e.reset(vcf_line);
		e.parse_basic_entry(true);

		if (e.get_N_alleles() > 2)
		{
			one_off_warning("\tIMPUTE: Only outputting biallelic loci.");
			continue;
		}

		// Exclude entries with missing data and/or unphased
		bool missing = false;
		for (ui=0; ui<N_indv; ui++)
		{
			if (include_indv[ui] == false)
				continue;

			if (include_genotype[s][ui] == false)
			{
				missing = true;
				break;
			}

			e.parse_genotype_entry(ui, true);
			e.get_indv_GENOTYPE_ids(ui, alleles);
			if ((alleles.first == -1) || (alleles.second == -1))
			{
				missing = true;
				break;
			}

			if (e.get_indv_PHASE(ui) != '|')
			{
				missing = true;
				break;
			}
		}
		if (missing == true)
			continue;

		if (e.get_ID() == ".")
		{
			legend << e.get_CHROM() << "-" << e.get_POS() << " " << e.get_POS() << " " << e.get_REF() << " " << e.get_ALT_allele(0) << endl;
		}
		else
			legend << e.get_ID() << " " << e.get_POS() << " " << e.get_REF() << " " << e.get_ALT_allele(0) << endl;

		bool first = true;
		for (ui=0; ui<N_indv; ui++)
		{
			if (include_indv[ui] == false)
				continue;

			e.parse_genotype_entry(ui, true);
			e.get_indv_GENOTYPE_ids(ui, alleles);
			if (first == true)
			{
				hap << alleles.first << " " << alleles.second;
				first = false;
			}
			else
				hap << " " << alleles.first << " " << alleles.second;
		}
		hap << endl;
	}

	hap.close();
	legend.close();
}

void vcf_file::output_LDhat_locs_file(const string &output_file_prefix, const string &chr, unsigned int &n_sites_out)
{
	string locs_file = output_file_prefix + ".ldhat.locs";
	ofstream locs(locs_file.c_str());
	if (!locs.is_open())
		error("Could not open LDhat locs Output File: " + locs_file, 2);

	int max_pos = -1;
	unsigned int n_sites=0;

	vcf_entry e(N_indv);
	string vcf_line;
	string chrom;
	for (unsigned int s=0; s<N_entries; s++)
	{
		if (include_entry[s] == false)
			continue;

		get_vcf_entry(s, vcf_line);
		e.reset(vcf_line);
		e.parse_basic_entry(true);

		if (e.get_N_alleles() != 2)
		{
			continue;
		}

		e.get_CHROM(chrom);
		if (chrom != chr)
			error("Mismatching chromosome in LDhat loci", 13);

		max_pos = max(e.get_POS(), max_pos);
		n_sites++;
	}

	locs << n_sites;
	locs.setf(ios::fixed,ios::floatfield);
	locs.precision(4);
	locs << "\t" << max_pos / 1000.0 << "\tL" << endl;
	for (unsigned int s=0; s<N_entries; s++)
	{
		if (include_entry[s] == false)
			continue;

		get_vcf_entry(s, vcf_line);
		e.reset(vcf_line);
		e.parse_basic_entry(true);

		if (e.get_N_alleles() != 2)
		{
			one_off_warning("\tLDhat: Only outputting biallelic loci.");
			continue;
		}

		locs << e.get_POS() / 1000.0 << endl;
	}
	locs.close();

	n_sites_out = n_sites;
}

void vcf_file::output_as_LDhat_phased(const string &output_file_prefix, const string &chr)
{
	if (has_genotypes == false)
		error("Require Genotypes in VCF file in order to output LDhat format.");

	printLOG("Outputting in phased LDhat format\n");
	if (chr == "")
		error("Require chromosome for LDhat output", 10);

	unsigned int n_sites;
	output_LDhat_locs_file(output_file_prefix, chr, n_sites);

	string sites_file = output_file_prefix + ".ldhat.sites";
	ofstream sites(sites_file.c_str());
	if (!sites.is_open())
		error("Could not open LDhat sites Output File: " + sites_file, 2);

	unsigned int n_indv = N_kept_individuals();
	pair<int, int> alleles;

	sites << n_indv*2 << "\t" << n_sites << "\t1" << endl;	// Note - this is incorrect for the X-chr.

	vector<ofstream *> tmp_files(2*N_indv);
	vector<string> tmp_filenames(2*N_indv);
	for (unsigned int ui=0; ui<N_indv; ui++)
	{
		if (include_indv[ui] == false)
			continue;
		string filename(tmpnam(NULL));
		ofstream *tmp_file = new ofstream(filename.c_str());
		if (!tmp_file->good())
			error("Could not open temp file.\n", 12);
		tmp_files[2*ui] = tmp_file;
		tmp_filenames[2*ui] = filename;

		string filename2(tmpnam(NULL));
		ofstream *tmp_file2 = new ofstream(filename2.c_str());
		if (!tmp_file2->good())
			error("\n\nCould not open temporary file.\n\n"
				"Most likely this is because the system is not allowing me to open enough temporary files.\n"
				"Try using ulimit -n <int> to increase the number of allowed open files.\n", 12);
		tmp_files[2*ui+1] = tmp_file2;
		tmp_filenames[2*ui+1] = filename2;
	}

	string vcf_line;
	vcf_entry e(N_indv);
	ofstream *tmp_file;

	for (unsigned int s=0; s<N_entries; s++)
	{
		if (include_entry[s] == false)
			continue;

		get_vcf_entry(s, vcf_line);
		e.reset(vcf_line);
		e.parse_basic_entry(true);

		if (e.get_N_alleles() != 2)
		{
			one_off_warning("\tLDhat: Only outputting biallelic loci.");
			continue;
		}

		for (unsigned int ui=0; ui<N_indv; ui++)
		{
			if (include_indv[ui] == false)
				continue;

			e.parse_genotype_entry(ui, true);
			e.get_indv_GENOTYPE_ids(ui, alleles);

			for (unsigned int k=0; k<2; k++)
			{
				tmp_file = tmp_files[(2*ui)+k];

				int geno;
				if (k == 0)
					geno = alleles.first;
				else
					geno = alleles.second;

				if ((geno != -1) && (include_genotype[s][ui]==true))
					(*tmp_file) << geno;
				else
					(*tmp_file) << "?";
			}
		}
	}

	string tmp_line;
	for (unsigned int ui=0; ui<N_indv; ui++)
	{
		if (include_indv[ui] == false)
			continue;

		for (unsigned int k=0; k<2; k++)
		{
			ofstream *tmp_file = tmp_files[2*ui+k];
			(*tmp_file) << endl;
			tmp_file->close();

			ifstream read_file(tmp_filenames[2*ui+k].c_str());
			if (!read_file.good())
				error("\n\nCould not open temporary file.\n\n"
				"Most likely this is because the system is not allowing me to open enough temporary files.\n"
				"Try using ulimit -n <int> to increase the number of allowed open files.\n", 12);
			getline(read_file, tmp_line);
			sites << ">" << indv[ui] << "-" << k << endl;
			sites << tmp_line << endl;
			read_file.close();
			remove(tmp_filenames[2*ui+k].c_str());
		}
	}

	sites.close();
}

void vcf_file::output_as_LDhat_unphased(const string &output_file_prefix, const string &chr)
{
	if (has_genotypes == false)
		error("Require Genotypes in VCF file in order to output LDhat format.");

	printLOG("Outputting in unphased LDhat format\n");
	if (chr == "")
		error("Require chromosome for LDhat output", 10);

	unsigned int n_sites;
	output_LDhat_locs_file(output_file_prefix, chr, n_sites);

	string sites_file = output_file_prefix + ".ldhat.sites";
	ofstream sites(sites_file.c_str());
	if (!sites.is_open())
		error("Could not open LDhat sites Output File: " + sites_file, 2);

	unsigned int n_indv = N_kept_individuals();
	pair<int, int> alleles;

	sites << n_indv << "\t" << n_sites << "\t2" << endl;

	vector<ofstream *> tmp_files(N_indv);
	vector<string> tmp_filenames(N_indv);
	for (unsigned int ui=0; ui<N_indv; ui++)
	{
		if (include_indv[ui] == false)
			continue;
		string filename(tmpnam(NULL));
		ofstream *tmp_file = new ofstream(filename.c_str());
		if (!tmp_file->good())
			error("\n\nCould not open temporary file.\n\n"
					"Most likely this is because the system is not allowing me to open enough temporary files.\n"
					"Try using ulimit -n <int> to increase the number of allowed open files.\n", 12);
		tmp_files[ui] = tmp_file;
		tmp_filenames[ui] = filename;
	}

	string vcf_line;
	vcf_entry e(N_indv);
	ofstream *tmp_file;
	for (unsigned int s=0; s<N_entries; s++)
	{
		if (include_entry[s] == false)
			continue;

		get_vcf_entry(s, vcf_line);
		e.reset(vcf_line);
		e.parse_basic_entry(true);

		if (e.get_N_alleles() != 2)
		{
			one_off_warning("\tLDhat: Only outputting biallelic loci.");
			continue;
		}

		for (unsigned int ui=0; ui<N_indv; ui++)
		{
			if (include_indv[ui] == false)
				continue;

			tmp_file = tmp_files[ui];

			if (include_genotype[s][ui] == false)
				(*tmp_file) << "?";
			else
			{
				e.parse_genotype_entry(ui, true);
				e.get_indv_GENOTYPE_ids(ui, alleles);

				switch (alleles.first)
				{
				case -1:
					(*tmp_file) << "?"; break;
				case 0:
					if (alleles.second == 0)
						(*tmp_file) << 0;
					else if (alleles.second == 1)
						(*tmp_file) << 2;
					else if ((alleles.second == -1) && (e.get_indv_PHASE(ui) == '|'))
						(*tmp_file) << 0;	// Haploid case
					else
						(*tmp_file) << '?';
					break;
				case 1:
					if (alleles.second == 0)
						(*tmp_file) << 2;
					else if (alleles.second == 1)
						(*tmp_file) << 1;
					else if ((alleles.second == -1) && (e.get_indv_PHASE(ui) == '|'))
						(*tmp_file) << 1;	// Haploid case
					else
						(*tmp_file) << '?';
					break;
				default:
					(*tmp_file) << '?';
				}
			}
		}
	}

	string tmp_line;
	for (unsigned int ui=0; ui<N_indv; ui++)
	{
		if (include_indv[ui] == false)
			continue;

		ofstream *tmp_file = tmp_files[ui];
		(*tmp_file) << endl;
		tmp_file->close();

		ifstream read_file(tmp_filenames[ui].c_str());
		if (!read_file.good())
			error("\n\nCould not open temporary file.\n\n"
				"Most likely this is because the system is not allowing me to open enough temporary files.\n"
				"Try using ulimit -n <int> to increase the number of allowed open files.\n", 12);
		getline(read_file, tmp_line);
		sites << ">" << indv[ui] << endl;
		sites << tmp_line << endl;
		read_file.close();
		remove(tmp_filenames[ui].c_str());
	}

	sites.close();
}

// Output INFO fields in tab-delimited format
void vcf_file::output_INFO_for_each_site(const string &output_file_prefix, const vector<string> &INFO_to_extract)
{
	if (INFO_to_extract.size() == 0)
		return;

	printLOG("Outputting INFO for each site\n");
	string output = output_file_prefix + ".INFO";
	ofstream out(output.c_str());
	if (!out.is_open())
		error("Could not open INFO Output File: " + output, 3);

	out << "CHROM\tPOS\tREF\tALT";
	for (unsigned int ui=0; ui<INFO_to_extract.size(); ui++)
		out << "\t" << INFO_to_extract[ui];
	out << endl;

	string vcf_line;
	vcf_entry e(N_indv);
	for (unsigned int s=0; s<N_entries; s++)
	{
		if (include_entry[s] == false)
			continue;

		get_vcf_entry(s, vcf_line);
		e.reset(vcf_line);
		e.parse_basic_entry(true, false, true);

		out << e.get_CHROM() << "\t" << e.get_POS() << "\t" << e.get_REF() << "\t" << e.get_ALT();

		for (unsigned int ui=0; ui<INFO_to_extract.size(); ui++)
		{
			out << "\t" << e.get_INFO_value(INFO_to_extract[ui]);
		}
		out << endl;
	}

	out.close();
}


// Output FORMAT information in tab-delimited format.
void vcf_file::output_FORMAT_information(const string &output_file_prefix, const string &FORMAT_id)
{
	if (FORMAT_id == "")
		return;

	if (has_genotypes == false)
		error("Require Genotypes in VCF file in order to output FORMAT information.");

	printLOG("Outputting FORMAT information for " + FORMAT_id + "\n");
	string output = output_file_prefix + "." + FORMAT_id + ".FORMAT";
	ofstream out(output.c_str());
	if (!out.is_open())
		error("Could not open FORMAT Output File: " + output, 7);

	out << "CHROM\tPOS";
	for (unsigned int ui=0; ui<N_indv; ui++)
	{
		if (include_indv[ui] == true)
			out << "\t" << indv[ui];
	}
	out << endl;

	string vcf_line, FORMAT_out;
	vcf_entry e(N_indv);
	for (unsigned int s=0; s<N_entries; s++)
	{
		if (include_entry[s] == false)
			continue;

		get_vcf_entry(s, vcf_line);
		e.reset(vcf_line);
		e.parse_basic_entry();
		e.parse_full_entry(true);

		if (e.FORMAT_id_exists(FORMAT_id) == false)
			continue;

		out << e.get_CHROM() << "\t" << e.get_POS();

		for (unsigned int ui=0; ui<N_indv; ui++)
		{
			if (include_indv[ui] == false)
				continue;

			e.read_indv_generic_entry(ui, FORMAT_id, FORMAT_out);
			out << "\t" << FORMAT_out;
		}
		out << endl;
	}
	out.close();
}

// Output genotype likelihoods from GL FORMAT tag, ready for input into BEAGLE
// using the Genotype likelihoods file format.
void vcf_file::output_BEAGLE_genotype_likelihoods(const string &output_file_prefix)
{
	if (has_genotypes == false)
		error("Require Genotypes in VCF file in order to output BEAGLE genotype likelihoods.");

	printLOG("Outputting in BEAGLE Genotype Likelihood format (bi-allelic SNPs with GL tags only)\n");

	string output = output_file_prefix + ".BEAGLE.GL";
	ofstream out(output.c_str());
	if (!out.is_open())
		error("Could not open BEAGLE GL Output File: " + output, 3);
	out << "marker\talleleA\talleleB";
	for (unsigned int ui=0; ui<N_indv; ui++)
	{
		if (include_indv[ui] == true)
			out << "\t" << indv[ui] << "\t" << indv[ui] << "\t" << indv[ui];
	}
	out << endl;

	string vcf_line, GL_entry, tmp_string;
	vcf_entry e(N_indv);
	double lk1, lk2, lk3;
	bool found_GL=false;
	istringstream ss;

	for (unsigned int s=0; s<N_entries; s++)
	{
		if (include_entry[s] == false)
			continue;

		get_vcf_entry(s, vcf_line);
		e.reset(vcf_line);
		e.parse_basic_entry(true);

		if (e.get_N_alleles() != 2)
		{
			one_off_warning("\tBEAGLE: Only outputting biallelic loci.");
			continue;
		}

		e.parse_full_entry(true);

		if (e.FORMAT_id_exists("GL") == false)
			continue;
		found_GL = true;

		out << e.get_CHROM() << ":" << e.get_POS() << "\t" << e.get_REF() << "\t" << e.get_ALT();

		for (unsigned int ui=0; ui<N_indv; ui++)
		{
			if (include_indv[ui] == false)
				continue;

			if (include_genotype[s][ui] == true)
			{
				e.read_indv_generic_entry(ui, "GL", GL_entry);
				ss.clear();
				ss.str(GL_entry);
				getline(ss, tmp_string, ',');
				lk1 = atof(tmp_string.c_str());
				getline(ss, tmp_string, ',');
				lk2 = atof(tmp_string.c_str());
				getline(ss, tmp_string);
				lk3 = atof(tmp_string.c_str());
				out << "\t" << pow(10,lk1) << "\t" << pow(10,lk2) << "\t" << pow(10,lk3);
			}
			else
			{
				out << "\t1\t1\t1";	// Mark as unknown
			}
		}
		out << endl;
	}

	if (found_GL == false)
		error("Require GL FORMAT tags in VCF file to output BEAGLE input.");
}

