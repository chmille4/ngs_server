/*
 * vcf_file_merge.cpp
 *
 *  Created on: Oct 30, 2009
 *      Author: Adam Auton
 *      ($Revision: 230 $)
 */

#include "vcf_file.h"

void vcf_file::return_site_union(vcf_file &file2, map<pair<string, int>, pair<int, int> > &CHROMPOS_to_filepos_pair)
{
	unsigned int s;
	int POS;
	string CHROM;
	string vcf_line;
	for (s=0; s<N_entries; s++)
	{
		if (include_entry[s] == true)
		{
			get_vcf_entry(s, vcf_line);
			vcf_entry e(N_indv, vcf_line);
			e.parse_basic_entry();

			CHROM = e.get_CHROM();
			POS = e.get_POS();

			CHROMPOS_to_filepos_pair[make_pair<string,int>(CHROM, POS)] = make_pair<int,int>(s, -1);
		}
	}
	for (s=0; s<file2.N_entries; s++)
	{
		if (file2.include_entry[s] == true)
		{
			file2.get_vcf_entry(s, vcf_line);
			vcf_entry e(file2.N_indv, vcf_line);
			e.parse_basic_entry();

			CHROM = e.get_CHROM();
			POS = e.get_POS();

			if (CHROMPOS_to_filepos_pair.find(make_pair<string,int>(CHROM, POS)) != CHROMPOS_to_filepos_pair.end())
			{
				CHROMPOS_to_filepos_pair[make_pair<string,int>(CHROM, POS)].second = s;
			}
			else
			{
				CHROMPOS_to_filepos_pair[make_pair<string,int>(CHROM, POS)] = make_pair<int,int>(-1, s);
			}
		}
	}
}


void vcf_file::return_indv_union(vcf_file &file2, map<string, pair< int, int> > &combined_individuals)
{
	for (unsigned int ui=0; ui<N_indv; ui++)
		if (include_indv[ui] == true)
			combined_individuals[indv[ui]] = make_pair<int,int>(ui, -1);

	for (unsigned int ui=0; ui<file2.N_indv; ui++)
		if (file2.include_indv[ui] == true)
		{
			if (combined_individuals.find(file2.indv[ui]) != combined_individuals.end())
				combined_individuals[file2.indv[ui]].second = ui;
			else
				combined_individuals[file2.indv[ui]] = make_pair<int,int>(-1, ui);
		}
}

void vcf_file::output_sites_in_files(const string &output_file_prefix, vcf_file &diff_vcf_file)
{
	printLOG("Comparing sites in VCF files...\n");
	map<pair<string, int>, pair<int, int> > CHROMPOS_to_filepos_pair;
	map<pair<string, int>, pair<int, int> >::iterator CHROMPOS_to_filepos_pair_it;
	return_site_union(diff_vcf_file, CHROMPOS_to_filepos_pair);

	string vcf_line;
	string CHROM;
	int POS;

	string output_file = output_file_prefix + ".diff.sites_in_files";
	ofstream sites_in_files(output_file.c_str());
	sites_in_files << "CHROM\tPOS\tIN_FILE\tREF\tALT1\tALT2" << endl;

	int s1, s2;
	int N_common_SNPs = 0, N_SNPs_file1_only=0, N_SNPs_file2_only=0;
	for (CHROMPOS_to_filepos_pair_it=CHROMPOS_to_filepos_pair.begin(); CHROMPOS_to_filepos_pair_it!=CHROMPOS_to_filepos_pair.end(); ++CHROMPOS_to_filepos_pair_it)
	{
		s1 = CHROMPOS_to_filepos_pair_it->second.first;
		s2 = CHROMPOS_to_filepos_pair_it->second.second;

		CHROM = CHROMPOS_to_filepos_pair_it->first.first;
		POS = CHROMPOS_to_filepos_pair_it->first.second;

		vcf_entry e1(N_indv);
		vcf_entry e2(diff_vcf_file.N_indv);

		// Read entries from file (if available)
		if (s1 != -1)
		{
			get_vcf_entry(s1, vcf_line);
			e1.reset(vcf_line);
		}

		if (s2 != -1)
		{
			diff_vcf_file.get_vcf_entry(s2, vcf_line);
			e2.reset(vcf_line);
		}

		e1.parse_basic_entry(true);
		e2.parse_basic_entry(true);

		// Set the reference to the non-missing entry (if available)
		string REF = e1.get_REF();
		string REF2 = e2.get_REF();
		if ((REF == "N") || (REF == "."))
			REF = REF2;
		if ((REF2 == "N") || (REF2 == "."))
			REF2 = REF;

		if ((REF != REF2) && (REF2 != "N") && (REF != "N") && (REF != ".") && (REF2 != "."))
			warning("Non-matching REF at " + CHROM + ":" + int2str(POS) + " " + REF + "/" + REF2 + ". Diff results may be unreliable.");

		sites_in_files << CHROM << "\t" << POS << "\t";
		if ((s1 != -1) && (s2 != -1))
		{
			N_common_SNPs++;
			sites_in_files << "B";
		}
		else if ((s1 != -1) && (s2 == -1))
		{
			N_SNPs_file1_only++;
			sites_in_files << "1";
		}
		else if ((s1 == -1) && (s2 != -1))
		{
			N_SNPs_file2_only++;
			sites_in_files << "2";
		}
		else
			error("SNP in neither file!?");

		sites_in_files << "\t" << REF << "\t" << e1.get_ALT() << "\t" << e2.get_ALT() << endl;
	}

	sites_in_files.close();

	printLOG("Found " + int2str(N_common_SNPs) + " SNPs common to both files.\n");
	printLOG("Found " + int2str(N_SNPs_file1_only) + " SNPs only in main file.\n");
	printLOG("Found " + int2str(N_SNPs_file2_only) + " SNPs only in second file.\n");
}

void vcf_file::output_indv_in_files(const string &output_file_prefix, vcf_file &diff_vcf_file)
{
	printLOG("Comparing individuals in VCF files...\n");

	string output_file = output_file_prefix + ".diff.indv_in_files";

	ofstream out(output_file.c_str());
	if (!out.is_open())
		error("Could not open Indv Differences File: " + output_file, 3);
	out << "INDV\tFILES" << endl;

	// Build a list of individuals contained in each file
	map<string, pair< int, int> > combined_individuals;
	map<string, pair< int, int> >::iterator combined_individuals_it;
	return_indv_union(diff_vcf_file, combined_individuals);

	unsigned int N_combined_indv = combined_individuals.size();
	unsigned int N[3]={0,0,0};
	for (combined_individuals_it=combined_individuals.begin(); combined_individuals_it!=combined_individuals.end(); ++combined_individuals_it)
	{
		if ((combined_individuals_it->second.first != -1) && (combined_individuals_it->second.second != -1))
		{
			N[0]++;
			out << combined_individuals_it->first << "\tB" << endl;
		}
		else if (combined_individuals_it->second.first != -1)
		{
			N[1]++;
			out << combined_individuals_it->first << "\t1" << endl;
		}
		else if (combined_individuals_it->second.second != -1)
		{
			N[2]++;
			out << combined_individuals_it->first << "\t2" << endl;
		}
		else
			error("Unhandled case");
	}
	out.close();

	printLOG("N_combined_individuals:\t" + int2str(N_combined_indv) + "\n");
	printLOG("N_individuals_common_to_both_files:\t" + int2str(N[0]) + "\n");
	printLOG("N_individuals_unique_to_file1:\t" + int2str(N[1]) + "\n");
	printLOG("N_individuals_unique_to_file2:\t" + int2str(N[2]) + "\n");
}

void vcf_file::output_discordance_by_indv(const string &output_file_prefix, vcf_file &diff_vcf_file)
{
	printLOG("Outputting Discordance By Individual...\n");
	map<pair<string, int>, pair<int, int> > CHROMPOS_to_filepos_pair;
	map<pair<string, int>, pair<int, int> >::iterator CHROMPOS_to_filepos_pair_it;
	return_site_union(diff_vcf_file, CHROMPOS_to_filepos_pair);

	map<string, pair< int, int> > combined_individuals;
	map<string, pair< int, int> >::iterator combined_individuals_it;
	return_indv_union(diff_vcf_file, combined_individuals);

	map<string, pair<int, int> > indv_sums;

	string vcf_line, CHROM;
	int POS;
	int s1, s2, indv1, indv2;

	for (CHROMPOS_to_filepos_pair_it=CHROMPOS_to_filepos_pair.begin(); CHROMPOS_to_filepos_pair_it != CHROMPOS_to_filepos_pair.end(); ++CHROMPOS_to_filepos_pair_it)
	{
		CHROM = CHROMPOS_to_filepos_pair_it->first.first;
		POS = CHROMPOS_to_filepos_pair_it->first.second;

		s1 = CHROMPOS_to_filepos_pair_it->second.first;
		s2 = CHROMPOS_to_filepos_pair_it->second.second;

		vcf_entry e1(N_indv);
		vcf_entry e2(diff_vcf_file.N_indv);

		// Read entries from file (if available)
		if (s1 != -1)
		{
			get_vcf_entry(s1, vcf_line);
			e1.reset(vcf_line);
		}

		if (s2 != -1)
		{
			diff_vcf_file.get_vcf_entry(s2, vcf_line);
			e2.reset(vcf_line);
		}

		e1.parse_basic_entry(true);
		e2.parse_basic_entry(true);

		// Set the reference to the non-missing entry (if available)
		string REF = e1.get_REF();
		string REF2 = e2.get_REF();
		if (REF == "N")
			REF = REF2;
		if (REF2 == "N")
			REF2 = REF;

		if (REF.size() != REF2.size())
		{
			warning("REF sequences at " + CHROM + ":" + int2str(POS) + " are not comparable. Skipping site");
			continue;
		}

		if ((REF != REF2) && (REF2 != "N") && (REF != "N"))
			warning("Non-matching REF " + CHROM + ":" + int2str(POS) + " " + REF + "/" + REF2);

		// Do the alternative alleles match?
		string ALT, ALT2;
		ALT = e1.get_ALT();
		ALT2 = e2.get_ALT();

		bool alleles_match = (ALT == ALT2) && (REF == REF2);
		e1.parse_full_entry(true);
		e1.parse_genotype_entries(true);

		e2.parse_full_entry(true);
		e2.parse_genotype_entries(true);

		pair<string, string> genotype1, genotype2;
		pair<int,int> geno_ids1, geno_ids2;
		pair<string, string> missing_genotype(".",".");
		pair<int, int> missing_id(-1,-1);

		for (combined_individuals_it=combined_individuals.begin(); combined_individuals_it!=combined_individuals.end(); ++combined_individuals_it)
		{
			indv1 = combined_individuals_it->second.first;
			indv2 = combined_individuals_it->second.second;

			if ((indv1 == -1) || (indv2 == -1))
				continue;	// Individual not found in one of the files

			if (alleles_match)
			{	// Alleles match, so can compare ids instead of strings
				e1.get_indv_GENOTYPE_ids(indv1, geno_ids1);
				e2.get_indv_GENOTYPE_ids(indv2, geno_ids2);

				if ((geno_ids1 != missing_id) && (geno_ids2 != missing_id))
				{
					indv_sums[combined_individuals_it->first].first++;
					if (((geno_ids1.first == geno_ids2.first) && (geno_ids1.second == geno_ids2.second)) ||
						((geno_ids1.first == geno_ids2.second) && (geno_ids1.second == geno_ids2.first)) )
					{	// Match
						// Don't do anything
					}
					else
					{	// Mismatch
						indv_sums[combined_individuals_it->first].second++;
					}
				}
				else if ((geno_ids1 == missing_id) && (geno_ids2 == missing_id))
				{	// Both missing
					// Don't do anything.
				}
				else if (geno_ids1 != missing_id)
				{	// Genotype 1 is not missing, genotype 2 is.
					// Don't do anything.
				}
				else if (geno_ids2 != missing_id)
				{	// Genotype 2 is not missing, genotype 1 is.
					// Don't do anything.
				}
				else
					error("Unknown condition");
			}
			else
			{	// Alleles don't match, so need to be more careful and compare strings
				e1.get_indv_GENOTYPE_strings(indv1, genotype1);
				e2.get_indv_GENOTYPE_strings(indv2, genotype2);

				if ((genotype1 != missing_genotype) && (genotype2 != missing_genotype))
				{	// No missing data
					indv_sums[combined_individuals_it->first].first++;
					if (((genotype1.first == genotype2.first) && (genotype1.second == genotype2.second)) ||
						((genotype1.first == genotype2.second) && (genotype1.second == genotype2.first)) )
					{	// Match
						// Don't do anything
					}
					else
					{	// Mismatch
						indv_sums[combined_individuals_it->first].second++;
					}
				}
				else if ((genotype1 == missing_genotype) && (genotype2 == missing_genotype))
				{	// Both missing
					// Don't do anything
				}
				else if (genotype1 != missing_genotype)
				{	// Genotype 1 is not missing, genotype 2 is.
					// Don't do anything
				}
				else if (genotype2 != missing_genotype)
				{	// Genotype 2 is not missing, genotype 1 is.
					// Don't do anything
				}
				else
					error("Unknown condition");
			}
		}
	}

	string output_file = output_file_prefix + ".diff.indv";
	ofstream out(output_file.c_str());
	if (!out.is_open())
		error("Could not open Sites Differences File: " + output_file, 3);
	out << "INDV\tN_COMMON_CALLED\tN_DISCORD\tDISCORDANCE" << endl;

	int N, N_discord;
	double discordance;
	for (combined_individuals_it=combined_individuals.begin(); combined_individuals_it!=combined_individuals.end(); ++combined_individuals_it)
	{
		out << combined_individuals_it->first;
		N = indv_sums[combined_individuals_it->first].first;
		N_discord = indv_sums[combined_individuals_it->first].second;
		discordance = N_discord / double(N);
		out << "\t" << N << "\t" << N_discord << "\t" << discordance << endl;
	}

	out.close();
}

void vcf_file::output_discordance_by_site(const string &output_file_prefix, vcf_file &diff_vcf_file)
{
	printLOG("Outputting Discordance By Site...\n");
	map<pair<string, int>, pair<int, int> > CHROMPOS_to_filepos_pair;
	map<pair<string, int>, pair<int, int> >::iterator CHROMPOS_to_filepos_pair_it;
	return_site_union(diff_vcf_file, CHROMPOS_to_filepos_pair);

	map<string, pair< int, int> > combined_individuals;
	map<string, pair< int, int> >::iterator combined_individuals_it;
	return_indv_union(diff_vcf_file, combined_individuals);

	string CHROM, vcf_line;
	int POS;
	int s1, s2, indv1, indv2;

	string output_file = output_file_prefix + ".diff.sites";
	ofstream diffsites(output_file.c_str());
	if (!diffsites.is_open())
		error("Could not open Sites Differences File: " + output_file, 3);
	//diffsites << "CHROM\tPOS\tFILES\tMATCHING_ALT\tN_COMMON_CALLED\tN_DISCORD\tDISCORDANCE\tN_FILE1_NONREF_GENOTYPES\tNON_REF_DISCORDANCE" << endl;
	diffsites << "CHROM\tPOS\tFILES\tMATCHING_ALLELES\tN_COMMON_CALLED\tN_DISCORD\tDISCORDANCE" << endl;

	for (CHROMPOS_to_filepos_pair_it=CHROMPOS_to_filepos_pair.begin(); CHROMPOS_to_filepos_pair_it != CHROMPOS_to_filepos_pair.end(); ++CHROMPOS_to_filepos_pair_it)
	{
		CHROM = CHROMPOS_to_filepos_pair_it->first.first;
		POS = CHROMPOS_to_filepos_pair_it->first.second;

		diffsites << CHROM << "\t" << POS;

		s1 = CHROMPOS_to_filepos_pair_it->second.first;
		s2 = CHROMPOS_to_filepos_pair_it->second.second;

		vcf_entry e1(N_indv);
		vcf_entry e2(diff_vcf_file.N_indv);

		bool data_in_both = true;
		// Read entries from file (if available)
		if (s1 != -1)
		{
			get_vcf_entry(s1, vcf_line);
			e1.reset(vcf_line);
		}
		else
			data_in_both = false;

		if (s2 != -1)
		{
			diff_vcf_file.get_vcf_entry(s2, vcf_line);
			e2.reset(vcf_line);
		}
		else
			data_in_both = false;

		if (data_in_both)
			diffsites << "\tB";
		else if ((s1 != -1) && (s2 == -1))
			diffsites << "\t1";
		else if ((s1 == -1) && (s2 != -1))
			diffsites << "\t2";
		else
			error("Unhandled condition");

		e1.parse_basic_entry(true);
		e2.parse_basic_entry(true);

		// Set the reference to the non-missing entry (if available)
		string REF = e1.get_REF();
		string REF2 = e2.get_REF();
		if (REF == "N")
			REF = REF2;
		if (REF2 == "N")
			REF2 = REF;

		if (REF.size() != REF2.size())
		{
			warning("REF sequences at " + CHROM + ":" + int2str(POS) + " are not comparable. Skipping site");
			continue;
		}

		if ((REF != REF2) && (REF2 != "N") && (REF != "N"))
			warning("Non-matching REF " + CHROM + ":" + int2str(POS) + " " + REF + "/" + REF2);

		// Do the alternative alleles match?
		string ALT, ALT2;
		ALT = e1.get_ALT();
		ALT2 = e2.get_ALT();

		bool alleles_match = ((ALT == ALT2) && (REF == REF2));
		diffsites << "\t" << alleles_match;

		e1.parse_full_entry(true);
		e1.parse_genotype_entries(true);

		e2.parse_full_entry(true);
		e2.parse_genotype_entries(true);

		pair<string, string> genotype1, genotype2;
		pair<int,int> geno_ids1, geno_ids2;
		pair<string, string> missing_genotype(".",".");
		pair<int, int> missing_id(-1,-1);

		unsigned int N_common_called=0;	// Number of genotypes called in both files
		unsigned int N_missing_1=0, N_missing_2=0;
		unsigned int N_discord=0;
		unsigned int N_concord_non_missing=0;

		for (combined_individuals_it=combined_individuals.begin(); combined_individuals_it!=combined_individuals.end(); ++combined_individuals_it)
		{
			indv1 = combined_individuals_it->second.first;
			indv2 = combined_individuals_it->second.second;

			if ((indv1 == -1) || (indv2 == -1))
				continue;	// Individual not found in one of the files

			if (alleles_match)
			{	// Alleles match, so can compare ids instead of strings
				e1.get_indv_GENOTYPE_ids(indv1, geno_ids1);
				e2.get_indv_GENOTYPE_ids(indv2, geno_ids2);

				if ((geno_ids1 != missing_id) && (geno_ids2 != missing_id))
				{
					N_common_called++;
					if (((geno_ids1.first == geno_ids2.first) && (geno_ids1.second == geno_ids2.second)) ||
						((geno_ids1.first == geno_ids2.second) && (geno_ids1.second == geno_ids2.first)) )
					{	// Match
						N_concord_non_missing++;
					}
					else
					{	// Mismatch
						N_discord++;
					}
				}
				else if ((geno_ids1 == missing_id) && (geno_ids2 == missing_id))
				{	// Both missing
					N_missing_1++; N_missing_2++;
				}
				else if (geno_ids1 != missing_id)
				{	// Genotype 1 is not missing, genotype 2 is.
					N_missing_2++;
				}
				else if (geno_ids2 != missing_id)
				{	// Genotype 2 is not missing, genotype 1 is.
					N_missing_1++;
				}
				else
					error("Unknown condition");
			}
			else
			{	// Alleles don't match, so need to be more careful and compare strings
				e1.get_indv_GENOTYPE_strings(indv1, genotype1);
				e2.get_indv_GENOTYPE_strings(indv2, genotype2);

				if ((genotype1 != missing_genotype) && (genotype2 != missing_genotype))
				{	// No missing data
					N_common_called++;
					if (((genotype1.first == genotype2.first) && (genotype1.second == genotype2.second)) ||
						((genotype1.first == genotype2.second) && (genotype1.second == genotype2.first)) )
					{	// Match
						N_concord_non_missing++;
					}
					else
					{	// Mismatch
						N_discord++;
					}
				}
				else if ((genotype1 == missing_genotype) && (genotype2 == missing_genotype))
				{	// Both missing
					N_missing_1++; N_missing_2++;
				}
				else if (genotype1 != missing_genotype)
				{	// Genotype 1 is not missing, genotype 2 is.
					N_missing_2++;
				}
				else if (genotype2 != missing_genotype)
				{	// Genotype 2 is not missing, genotype 1 is.
					N_missing_1++;
				}
				else
					error("Unknown condition");
			}
		}
		double discordance = N_discord / double(N_common_called);
		diffsites << "\t" << N_common_called << "\t" << N_discord << "\t" << discordance;
		diffsites << endl;
	}
	diffsites.close();
}

void vcf_file::output_discordance_matrix(const string &output_file_prefix, vcf_file &diff_vcf_file)
{
	printLOG("Outputting Discordance Matrix\n\tFor bi-allelic loci, called in both files, with matching alleles only...\n");
	map<pair<string, int>, pair<int, int> > CHROMPOS_to_filepos_pair;
	map<pair<string, int>, pair<int, int> >::iterator CHROMPOS_to_filepos_pair_it;
	return_site_union(diff_vcf_file, CHROMPOS_to_filepos_pair);

	map<string, pair< int, int> > combined_individuals;
	map<string, pair< int, int> >::iterator combined_individuals_it;
	return_indv_union(diff_vcf_file, combined_individuals);

	string vcf_line;
	int s1, s2, indv1, indv2;

	vector<vector<int> > discordance_matrix(4, vector<int>(4, 0));

	for (CHROMPOS_to_filepos_pair_it=CHROMPOS_to_filepos_pair.begin(); CHROMPOS_to_filepos_pair_it != CHROMPOS_to_filepos_pair.end(); ++CHROMPOS_to_filepos_pair_it)
	{
		s1 = CHROMPOS_to_filepos_pair_it->second.first;
		s2 = CHROMPOS_to_filepos_pair_it->second.second;

		vcf_entry e1(N_indv);
		vcf_entry e2(diff_vcf_file.N_indv);

		// Read entries from file (if available)
		if (s1 != -1)
		{
			get_vcf_entry(s1, vcf_line);
			e1.reset(vcf_line);
		}

		if (s2 != -1)
		{
			diff_vcf_file.get_vcf_entry(s2, vcf_line);
			e2.reset(vcf_line);
		}

		e1.parse_basic_entry(true);
		e2.parse_basic_entry(true);

		if ((e1.get_N_alleles() != 2) || (e2.get_N_alleles() != 2))
			continue;

		// Set the reference to the non-missing entry (if available)
		string REF = e1.get_REF();
		string REF2 = e2.get_REF();
		if (REF == "N")
			REF = REF2;
		if (REF2 == "N")
			REF2 = REF;

		if (REF.size() != REF2.size())
			continue;

		if ((REF != REF2) && (REF2 != "N") && (REF != "N"))
			continue;

		// Do the alternative alleles match?
		string ALT, ALT2;
		ALT = e1.get_ALT();
		ALT2 = e2.get_ALT();

		bool alleles_match = (ALT == ALT2) && (REF == REF2);
		if (alleles_match == false)
			continue;

		e1.parse_full_entry(true);
		e1.parse_genotype_entries(true);

		e2.parse_full_entry(true);
		e2.parse_genotype_entries(true);

		pair<int,int> geno_ids1, geno_ids2;
		int N1, N2;

		for (combined_individuals_it=combined_individuals.begin(); combined_individuals_it!=combined_individuals.end(); ++combined_individuals_it)
		{
			indv1 = combined_individuals_it->second.first;
			indv2 = combined_individuals_it->second.second;

			if ((indv1 == -1) || (indv2 == -1))
				continue;	// Individual not found in one of the files

			// Alleles match, so can compare ids instead of strings
			e1.get_indv_GENOTYPE_ids(indv1, geno_ids1);
			e2.get_indv_GENOTYPE_ids(indv2, geno_ids2);

			if (((geno_ids1.first != -1) && (geno_ids1.second == -1)) ||
				((geno_ids2.first != -1) && (geno_ids2.second == -1)))
			{	// Haploid
				one_off_warning("***Warning: Haploid chromosomes not counted!***");
				continue;
			}

			N1 = geno_ids1.first + geno_ids1.second;
			N2 = geno_ids2.first + geno_ids2.second;

			if ((N1 == -1) || (N1 < -2) || (N1 > 2))
				error("Unhandled case");
			if ((N2 == -1) || (N2 < -2) || (N2 > 2))
				error("Unhandled case");

			if (N1 == -2)
				N1 = 3;

			if (N2 == -2)
				N2 = 3;

			discordance_matrix[N1][N2]++;
		}
	}

	string output_file = output_file_prefix + ".diff.discordance_matrix";
	ofstream out(output_file.c_str());
	if (!out.is_open())
		error("Could not open Discordance Matrix File: " + output_file, 3);

	out << "-\tN_0/0_file1\tN_0/1_file1\tN_1/1_file1\tN_./._file1" << endl;
	out << "N_0/0_file2\t" << discordance_matrix[0][0] << "\t" << discordance_matrix[1][0] << "\t" << discordance_matrix[2][0] << "\t" << discordance_matrix[3][0] << endl;
	out << "N_0/1_file2\t" << discordance_matrix[0][1] << "\t" << discordance_matrix[1][1] << "\t" << discordance_matrix[2][1] << "\t" << discordance_matrix[3][1] << endl;
	out << "N_1/1_file2\t" << discordance_matrix[0][2] << "\t" << discordance_matrix[1][2] << "\t" << discordance_matrix[2][2] << "\t" << discordance_matrix[3][2] << endl;
	out << "N_./._file2\t" << discordance_matrix[0][3] << "\t" << discordance_matrix[1][3] << "\t" << discordance_matrix[2][3] << "\t" << discordance_matrix[3][3] << endl;
	out.close();
}

void vcf_file::output_switch_error(const string &output_file_prefix, vcf_file &diff_vcf_file)
{
	printLOG("Outputting Phase Switch Errors...\n");
	map<pair<string, int>, pair<int, int> > CHROMPOS_to_filepos_pair;
	map<pair<string, int>, pair<int, int> >::iterator CHROMPOS_to_filepos_pair_it;
	return_site_union(diff_vcf_file, CHROMPOS_to_filepos_pair);

	map<string, pair< int, int> > combined_individuals;
	map<string, pair< int, int> >::iterator combined_individuals_it;
	return_indv_union(diff_vcf_file, combined_individuals);

	string CHROM, vcf_line;
	int POS;
	int s1, s2, indv1, indv2;

	string output_file = output_file_prefix + ".diff.switch";
	ofstream switcherror(output_file.c_str());
	if (!switcherror.is_open())
		error("Could not open Switch Error file: " + output_file, 4);
	switcherror << "CHROM\tPOS\tINDV" << endl;

	unsigned int N_combined_indv = combined_individuals.size();
	vector<int> N_phased_het_sites(N_combined_indv, 0);
	vector<int> N_switch_errors(N_combined_indv, 0);

	pair<string, string> missing_genotype(".",".");
	vector<pair<string, string> > prev_geno_file1(N_combined_indv, missing_genotype);
	vector<pair<string, string> > prev_geno_file2(N_combined_indv, missing_genotype);
	pair<string, string> file1_hap1, file1_hap2, file2_hap1;

	for (CHROMPOS_to_filepos_pair_it=CHROMPOS_to_filepos_pair.begin(); CHROMPOS_to_filepos_pair_it != CHROMPOS_to_filepos_pair.end(); ++CHROMPOS_to_filepos_pair_it)
	{
		CHROM = CHROMPOS_to_filepos_pair_it->first.first;
		POS = CHROMPOS_to_filepos_pair_it->first.second;

		s1 = CHROMPOS_to_filepos_pair_it->second.first;
		s2 = CHROMPOS_to_filepos_pair_it->second.second;

		vcf_entry e1(N_indv);
		vcf_entry e2(diff_vcf_file.N_indv);

		// Read entries from file (if available)
		if (s1 != -1)
		{
			get_vcf_entry(s1, vcf_line);
			e1.reset(vcf_line);
		}

		if (s2 != -1)
		{
			diff_vcf_file.get_vcf_entry(s2, vcf_line);
			e2.reset(vcf_line);
		}

		e1.parse_basic_entry(true);
		e2.parse_basic_entry(true);

		e1.parse_full_entry(true);
		e1.parse_genotype_entries(true);

		e2.parse_full_entry(true);
		e2.parse_genotype_entries(true);

		pair<string, string> genotype1, genotype2;
		pair<string, string> missing_genotype(".",".");

		unsigned int N_common_called=0;	// Number of genotypes called in both files
		unsigned int indv_count=0;

		// Bug fix applied (#3354189) - July 5th 2011
		for (combined_individuals_it=combined_individuals.begin();
				combined_individuals_it!=combined_individuals.end();
				++combined_individuals_it, indv_count++)
		{
			indv1 = combined_individuals_it->second.first;
			indv2 = combined_individuals_it->second.second;

			if ((indv1 == -1) || (indv2 == -1))
				continue;	// Individual not found in one of the files

			e1.get_indv_GENOTYPE_strings(indv1, genotype1);
			e2.get_indv_GENOTYPE_strings(indv2, genotype2);

			if ((genotype1 != missing_genotype) && (genotype2 != missing_genotype))
			{	// No missing data
				N_common_called++;
				if (((genotype1.first == genotype2.first) && (genotype1.second == genotype2.second)) ||
					((genotype1.first == genotype2.second) && (genotype1.second == genotype2.first)) )
				{	// Have a matching genotypes in files 1 and 2
					if (genotype1.first != genotype1.second)
					{	// It's a heterozgote
						char phase1, phase2;
						phase1 = e1.get_indv_PHASE(indv1);
						phase2 = e2.get_indv_PHASE(indv2);
						if ((phase1 == '|') && (phase2 == '|'))
						{	// Calculate Phasing error (switch error)
							N_phased_het_sites[indv_count]++;
							file1_hap1 = make_pair<string,string>(prev_geno_file1[indv_count].first, genotype1.first);
							file1_hap2 = make_pair<string,string>(prev_geno_file1[indv_count].second, genotype1.second);
							file2_hap1 = make_pair<string,string>(prev_geno_file2[indv_count].first, genotype2.first);

							if ((file2_hap1 != file1_hap1) && (file2_hap1 != file1_hap2))
							{	// Must be a switch error
								string indv_id;
								N_switch_errors[indv_count]++;
								if (indv1 != -1)
									indv_id = indv[indv1];
								else
									indv_id = diff_vcf_file.indv[indv2];
								switcherror << CHROM << "\t" << POS << "\t" << indv_id << endl;
							}
							prev_geno_file1[indv_count] = genotype1;
							prev_geno_file2[indv_count] = genotype2;
						}
					}
				}
			}
		}
	}
	switcherror.close();

	output_file = output_file_prefix + ".diff.indv.switch";
	ofstream idiscord(output_file.c_str());
	if (!idiscord.is_open())
		error("Could not open Individual Discordance File: " + output_file, 3);

	idiscord << "INDV\tN_COMMON_PHASED_HET\tN_SWITCH\tSWITCH" << endl;
	unsigned int indv_count=0;
	double switch_error;
	string indv_id;
	for (combined_individuals_it=combined_individuals.begin(); combined_individuals_it!=combined_individuals.end(); ++combined_individuals_it)
	{
		indv1 = combined_individuals_it->second.first;
		indv2 = combined_individuals_it->second.second;

		if (indv1 != -1)
			indv_id = indv[indv1];
		else
			indv_id = diff_vcf_file.indv[indv2];

		if (N_phased_het_sites[indv_count] > 0)
			switch_error = double(N_switch_errors[indv_count]) / N_phased_het_sites[indv_count];
		else
			switch_error = 0;
		idiscord << indv_id << "\t" << N_phased_het_sites[indv_count] << "\t" << N_switch_errors[indv_count] << "\t" << switch_error << endl;

		indv_count++;
	}
	idiscord.close();
}

/*
void vcf_file::output_concensus_statistics(const string &output_file_prefix, vcf_file &diff_vcf_file)
{
	printLOG("Outputting Consensus Statistics... \n");
	unsigned int ui;

	string output_file;
	string vcf_line;

	// Build a list of individuals contained in each file
	map<string, pair< int, int> > combined_individuals;
	map<string, pair< int, int> >::iterator combined_individuals_it;
	for (ui=0; ui<N_indv; ui++)
		if (include_indv[ui] == true)
			combined_individuals[indv[ui]] = make_pair<int,int>(ui, -1);

	for (ui=0; ui<diff_vcf_file.N_indv; ui++)
		if (diff_vcf_file.include_indv[ui] == true)
		{
			if (combined_individuals.find(diff_vcf_file.indv[ui]) != combined_individuals.end())
				combined_individuals[diff_vcf_file.indv[ui]].second = ui;
			else
				combined_individuals[diff_vcf_file.indv[ui]] = make_pair<int,int>(-1, ui);
		}

	unsigned int N_combined_indv = combined_individuals.size();
	unsigned int N[3]={0,0,0};
	for (combined_individuals_it=combined_individuals.begin(); combined_individuals_it!=combined_individuals.end(); ++combined_individuals_it)
	{
		if ((combined_individuals_it->second.first != -1) && (combined_individuals_it->second.second != -1))
			N[0]++;
		else if (combined_individuals_it->second.first != -1)
			N[1]++;
		else
			N[2]++;
	}

	vector<int> indv_N_discord(N_combined_indv, 0);
	vector<int> indv_N_called_sites(N_combined_indv, 0);

	printLOG("N_combined_individuals:\t" + int2str(N_combined_indv) + "\n");
	printLOG("N_individuals_common_to_both_files:\t" + int2str(N[0]) + "\n");
	printLOG("N_individuals_unique_to_file1:\t" + int2str(N[1]) + "\n");
	printLOG("N_individuals_unique_to_file2:\t" + int2str(N[2]) + "\n");

	// Build a table of included entries in both files
	map<pair<string, int>, pair<int, int> > CHROMPOS_to_filepos_pair;
	map<pair<string, int>, pair<int, int> >::iterator CHROMPOS_to_filepos_pair_it;
	string CHROM;
	int POS;
	return_site_union(diff_vcf_file, CHROMPOS_to_filepos_pair);

	output_file = output_file_prefix + ".diff.sites_in_files";
	ofstream sites_in_files(output_file.c_str());
	sites_in_files << "CHROM\tPOS\tIN_FILE" << endl;

	int s1, s2;
	int N_common_SNPs = 0, N_SNPs_file1_only=0, N_SNPs_file2_only=0;
	for (CHROMPOS_to_filepos_pair_it=CHROMPOS_to_filepos_pair.begin(); CHROMPOS_to_filepos_pair_it!=CHROMPOS_to_filepos_pair.end(); ++CHROMPOS_to_filepos_pair_it)
	{
		s1 = CHROMPOS_to_filepos_pair_it->second.first;
		s2 = CHROMPOS_to_filepos_pair_it->second.second;

		sites_in_files << CHROMPOS_to_filepos_pair_it->first.first << "\t" << CHROMPOS_to_filepos_pair_it->first.second << "\t";
		if ((s1 != -1) && (s2 != -1))
		{
			N_common_SNPs++;
			sites_in_files << "B" << endl;
		}
		else if ((s1 != -1) && (s2 == -1))
		{
			N_SNPs_file1_only++;
			sites_in_files << "1" << endl;
		}
		else if ((s1 == -1) && (s2 != -1))
		{
			N_SNPs_file2_only++;
			sites_in_files << "2" << endl;
		}
		else
			error("SNP in neither file!?");
	}

	sites_in_files.close();

	printLOG("Found " + int2str(N_common_SNPs) + " SNPs common to both files.\n");
	printLOG("Found " + int2str(N_SNPs_file1_only) + " SNPs only in main file.\n");
	printLOG("Found " + int2str(N_SNPs_file2_only) + " SNPs only in second file.\n");

	output_file = output_file_prefix + ".diff.sites";
	ofstream diffsites(output_file.c_str());
	if (!diffsites.is_open())
		error("Could not open Sites Differences File: " + output_file, 3);
	diffsites << "CHROM\tPOS\tFILES\tMATCHING_ALT\tN_COMMON_CALLED\tN_DISCORD\tDISCORDANCE\tN_FILE1_NONREF_GENOTYPES\tNON_REF_DISCORDANCE" << endl;

	output_file = output_file_prefix + ".diff.switch";
	ofstream switcherror(output_file.c_str());
	if (!switcherror.is_open())
		error("Could not open Switch Error file: " + output_file, 4);
	switcherror << "CHROM\tPOS\tINDV" << endl;

	// Now try and merge the entries.
	unsigned int N_common_genotypes = 0;
	unsigned int N_common_discordant_genotypes = 0;
	unsigned int N_sites_with_mismatching_ALT = 0;
	unsigned int N_non_ref_genotypes = 0;
	unsigned int N_discordant_non_ref_genotypes = 0;

	pair<string, string> genotype1, genotype2;
	pair<int,int> geno_ids1, geno_ids2;
	pair<string, string> missing_genotype(".",".");
	pair<int, int> missing_HQUAL(0,0);

	pair<int, int> homo_ref(0, 0);
	vector<pair<string, string> > prev_geno_file1(N_combined_indv, missing_genotype);
	vector<pair<string, string> > prev_geno_file2(N_combined_indv, missing_genotype);
	pair<string, string> file1_hap1, file1_hap2, file2_hap1;

	vector<int> N_phased_het_sites(N_combined_indv, 0);
	vector<int> N_switch_errors(N_combined_indv, 0);
	vector<pair<int,int> > indv_depth_at_common_sites(N_combined_indv, make_pair(0,0));
	vector<pair<int,int> > indv_count_at_common_sites(N_combined_indv, make_pair(0,0));

	vector<vector<int> > genotype_concord_matrix(4, vector<int>(4, 0));

	for (CHROMPOS_to_filepos_pair_it=CHROMPOS_to_filepos_pair.begin(); CHROMPOS_to_filepos_pair_it != CHROMPOS_to_filepos_pair.end(); ++CHROMPOS_to_filepos_pair_it)
	{
		CHROM = CHROMPOS_to_filepos_pair_it->first.first;
		POS = CHROMPOS_to_filepos_pair_it->first.second;

		s1 = CHROMPOS_to_filepos_pair_it->second.first;
		s2 = CHROMPOS_to_filepos_pair_it->second.second;

		vcf_entry e1(N_indv);
		vcf_entry e2(diff_vcf_file.N_indv);

		bool data_in_both = true;
		// Read entries from file (if available)
		if (s1 != -1)
		{
			get_vcf_entry(s1, vcf_line);
			e1.reset(vcf_line);
		}
		else
			data_in_both = false;

		if (s2 != -1)
		{
			diff_vcf_file.get_vcf_entry(s2, vcf_line);
			e2.reset(vcf_line);
		}
		else
			data_in_both = false;

		e1.parse_basic_entry(true, true, true);
		e2.parse_basic_entry(true, true, true);

		// Set the reference to the non-missing entry (if available)
		string REF = e1.get_REF();
		string REF2 = e2.get_REF();
		if (REF == "N")
			REF = REF2;
		if (REF2 == "N")
			REF2 = REF;

		if (REF.size() != REF2.size())
		{
			warning("REF sequences at " + CHROM + ":" + int2str(POS) + " are not comparable. Skipping site");
			continue;
		}

		if ((REF != REF2) && (REF2 != "N") && (REF != "N"))
			warning("Non-matching REF " + CHROM + ":" + int2str(POS) + " " + REF + "/" + REF2);

		// Do the alternative alleles match?
		set<string> ALT1, ALT2;
		for (ui=0; ui<(e1.get_N_alleles()-1); ui++)
			ALT1.insert(e1.get_ALT_allele(ui));

		for (ui=0; ui<(e2.get_N_alleles()-1); ui++)
			ALT2.insert(e2.get_ALT_allele(ui));

		bool matching_ALT=true;
		if ((data_in_both) && (ALT1 != ALT2) && (ALT1.size() > 0) && (ALT2.size() > 0))
		{
			N_sites_with_mismatching_ALT++;
			matching_ALT = false;
		}

		if (data_in_both)
		{
			diffsites << CHROM << "\t" << POS << "\t";
			diffsites << "B\t" << matching_ALT << "\t";
		}
		else
		{
			continue;
		}

		if (s1 != -1)
		{
			e1.parse_full_entry(true);
			e1.parse_genotype_entries(true, true, true, true, true);
		}

		if (s2 != -1)
		{
			e2.parse_full_entry(true);
			e2.parse_genotype_entries(true, true, true, true, true);
		}

		// Now merge the genotypes.
		unsigned int indv_count=0;
		int indv1, indv2;
		unsigned int N_discordant_site_counter=0;
		unsigned int N_indvs_with_data=0;
		unsigned int site_N_non_ref_genotypes=0;
		unsigned int site_N_discordant_non_ref_genotypes = 0;
		int depth;
		for (combined_individuals_it=combined_individuals.begin(); combined_individuals_it!=combined_individuals.end(); ++combined_individuals_it)
		{
			indv1 = combined_individuals_it->second.first;
			indv2 = combined_individuals_it->second.second;

			if ((indv1 == -1) && (indv2 == -1))
			{	// Genotype is completely missing... should never happen
				error("Missing genotype!?", 83);
			}
			else if ((indv1 == -1) && (indv2 != -1))
			{	// Data is missing from first file, so just use second file.

			}
			else if ((indv1 != -1) && (indv2 == -1))
			{	// Data is missing from second file, so just use first file.

			}
			else
			{	// Data from both files, so figure out what to do
				bool non_ref_genotype = false;
				if (data_in_both)
				{
					e1.get_indv_GENOTYPE_strings(indv1, genotype1);
					e2.get_indv_GENOTYPE_strings(indv2, genotype2);
					e1.get_indv_GENOTYPE_ids(indv1, geno_ids1);
					e2.get_indv_GENOTYPE_ids(indv2, geno_ids2);

					N_common_genotypes++;
					N_indvs_with_data++;

					if (geno_ids1 != homo_ref)
					{	// First file is not a hom ref
						N_non_ref_genotypes++;
						site_N_non_ref_genotypes++;
						non_ref_genotype = true;
					}

					depth = e1.get_indv_DEPTH(indv1);
					if (depth >= 0)
					{
						indv_depth_at_common_sites[indv_count].first += depth;
						indv_count_at_common_sites[indv_count].first++;
					}
					depth = e2.get_indv_DEPTH(indv2);
					if (depth >= 0)
					{
						indv_depth_at_common_sites[indv_count].second += depth;
						indv_count_at_common_sites[indv_count].second++;
					}
				}

				if ((genotype1 == missing_genotype) && (genotype2 == missing_genotype))
				{
					genotype_concord_matrix[3][3]++;
				}

				if ((genotype1 == missing_genotype) && (genotype2 != missing_genotype))
				{	// Missing data, Favour second file
					if (matching_ALT && (ALT2.size() <= 1))
					{
						unsigned int idx2 = geno_ids2.first + geno_ids2.second;
						genotype_concord_matrix[3][idx2]++;
					}
				}

				if ((genotype2 == missing_genotype) && (genotype1 != missing_genotype))
				{	// Favour first file
					if (matching_ALT && (ALT1.size() <= 1))
					{
						unsigned int idx1 = geno_ids1.first + geno_ids1.second;
						genotype_concord_matrix[idx1][3]++;
					}
				}

				if ((genotype1 != missing_genotype) && (genotype2 != missing_genotype))
				{
					if (data_in_both)
					{
						if (matching_ALT && (ALT1.size() <= 1) && (ALT2.size() <= 1))
						{
							unsigned int idx1 = geno_ids1.first + geno_ids1.second;
							unsigned int idx2 = geno_ids2.first + geno_ids2.second;
							genotype_concord_matrix[idx1][idx2]++;
						}

						indv_N_called_sites[indv_count]++;
						if (!vcf_entry::genotypes_equal(genotype1, genotype2))
						{
							N_common_discordant_genotypes++;
							N_discordant_site_counter++;
							indv_N_discord[indv_count]++;

							if (non_ref_genotype)
							{
								N_discordant_non_ref_genotypes++;
								site_N_discordant_non_ref_genotypes++;
							}
						}
						else
						{	// Have a matching genotype in files 1 and 2
							if (geno_ids1.first != geno_ids1.second)
							{	// It's a heterozgote
								char phase1, phase2;
								phase1 = e1.get_indv_PHASE(indv1);
								phase2 = e2.get_indv_PHASE(indv2);
								if ((phase1 == '|') && (phase2 == '|'))
								{	// Calculate Phasing error (switch error)
									N_phased_het_sites[indv_count]++;
									file1_hap1 = make_pair<string,string>(prev_geno_file1[indv_count].first, genotype1.first);
									file1_hap2 = make_pair<string,string>(prev_geno_file1[indv_count].second, genotype1.second);
									file2_hap1 = make_pair<string,string>(prev_geno_file2[indv_count].first, genotype2.first);

									if ((file2_hap1 != file1_hap1) && (file2_hap1 != file1_hap2))
									{	// Must be a switch error
										string indv_id;
										N_switch_errors[indv_count]++;
										if (indv1 != -1)
											indv_id = indv[indv1];
										else
											indv_id = diff_vcf_file.indv[indv2];
										switcherror << CHROM << "\t" << POS << "\t" << indv_id << endl;
									}
									prev_geno_file1[indv_count] = genotype1;
									prev_geno_file2[indv_count] = genotype2;
								}
							}
						}
					}
				}
			}

			indv_count++;
		}
		double discordance = 0.0;
		if (N_indvs_with_data > 0)
			discordance = double(N_discordant_site_counter) / N_indvs_with_data;
		double non_ref_discordance = 0.0;
		if (site_N_non_ref_genotypes > 0)
			non_ref_discordance = double(site_N_discordant_non_ref_genotypes) / site_N_non_ref_genotypes;
		diffsites << N_indvs_with_data << "\t" << N_discordant_site_counter << "\t" << discordance;
		diffsites << "\t" << site_N_non_ref_genotypes << "\t" << non_ref_discordance;
		diffsites << endl;
	}

	output_file = output_file_prefix + ".diff.4x4";
	ofstream four_by_four(output_file.c_str());
	if (!four_by_four.is_open())
		error("Could not open 3x3 File: " + output_file, 3);

	four_by_four << "-\tN00_file1\tN01_file1\tN11_file1\tN.._file1" << endl;

	four_by_four << "N00_file2\t" << genotype_concord_matrix[0][0] << "\t" << genotype_concord_matrix[1][0] << "\t" << genotype_concord_matrix[2][0] << "\t" << genotype_concord_matrix[3][0] << endl;
	four_by_four << "N01_file2\t" << genotype_concord_matrix[0][1] << "\t" << genotype_concord_matrix[1][1] << "\t" << genotype_concord_matrix[2][1] << "\t" << genotype_concord_matrix[3][1] << endl;
	four_by_four << "N11_file2\t" << genotype_concord_matrix[0][2] << "\t" << genotype_concord_matrix[1][2] << "\t" << genotype_concord_matrix[2][2] << "\t" << genotype_concord_matrix[3][2] << endl;
	four_by_four << "N.._file2\t" << genotype_concord_matrix[0][3] << "\t" << genotype_concord_matrix[1][3] << "\t" << genotype_concord_matrix[2][3] << "\t" << genotype_concord_matrix[3][3] << endl;
	four_by_four.close();


	output_file = output_file_prefix + ".diff.indv.discord";
	ofstream idiscord(output_file.c_str());
	if (!idiscord.is_open())
		error("Could not open Individual Discordance File: " + output_file, 3);

	idiscord << "INDV\tMEAN_DP_1\tMEAN_DP_2\tN_COMMON_CALLED\tN_DISCORD\tDISCORD\tN_COMMON_PHASED_HET\tN_SWITCH\tSWITCH" << endl;
	unsigned int indv_count=0;
	double discordance, switch_error;
	int indv1, indv2;
	string indv_id;
	for (combined_individuals_it=combined_individuals.begin(); combined_individuals_it!=combined_individuals.end(); ++combined_individuals_it)
	{
		indv1 = combined_individuals_it->second.first;
		indv2 = combined_individuals_it->second.second;

		if (indv1 != -1)
			indv_id = indv[indv1];
		else
			indv_id = diff_vcf_file.indv[indv2];

		if (indv_N_called_sites[indv_count] > 0)
			discordance = double(indv_N_discord[indv_count]) / indv_N_called_sites[indv_count];
		else
			discordance = 0.0;
		idiscord << indv_id;

		double mean_depth1 = 0, mean_depth2=0;
		if (indv_count_at_common_sites[indv_count].first > 0)
		{
			mean_depth1 = double(indv_depth_at_common_sites[indv_count].first) / indv_count_at_common_sites[indv_count].first;
		}

		if (indv_count_at_common_sites[indv_count].second > 0)
		{
			mean_depth2 = double(indv_depth_at_common_sites[indv_count].second) / indv_count_at_common_sites[indv_count].second;
		}
		idiscord << "\t" << mean_depth1 << "\t" << mean_depth2;

		idiscord << "\t" << indv_N_called_sites[indv_count] << "\t" << indv_N_discord[indv_count] << "\t" << discordance;
		if (N_phased_het_sites[indv_count] > 0)
			switch_error = double(N_switch_errors[indv_count]) / N_phased_het_sites[indv_count];
		else
			switch_error = 0;
		idiscord << "\t" << N_phased_het_sites[indv_count] << "\t" << N_switch_errors[indv_count] << "\t" << switch_error << endl;

		indv_count++;
	}
	idiscord.close();

	printLOG("Found " + int2str(N_sites_with_mismatching_ALT) + " sites with mismatching ALT alleles.\n");

	printLOG("Found " + int2str(N_non_ref_genotypes) + " non-reference genotypes called in both files.\n");
	printLOG("Found " + int2str(N_discordant_non_ref_genotypes) + " discordant non-reference genotypes.\n");
	double concordance = 1.0 - (double(N_discordant_non_ref_genotypes)) / N_non_ref_genotypes;
	printLOG("Concordance rate: " + dbl2str_fixed(concordance * 100,2) + "%\n");

	printLOG("Found " + int2str(N_common_genotypes) + " genotypes called in both files.\n");
	printLOG("Found " + int2str(N_common_discordant_genotypes) + " discordant genotypes.\n");
	concordance = 1.0 - (double(N_common_discordant_genotypes)) / N_common_genotypes;
	printLOG("Overall Concordance rate: " + dbl2str_fixed(concordance * 100,2) + "%\n");

	diffsites.close();
	switcherror.close();
	printLOG("Done\n");
}
*/
