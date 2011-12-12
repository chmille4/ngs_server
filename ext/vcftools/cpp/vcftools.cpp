/*
 * vcftools.cpp
 *
 *  Created on: Aug 19, 2009
 *      Author: Adam Auton
 *      ($Revision: 249 $)
 */
#include "vcftools.h"

ofstream LOG;

int main(int argc, char *argv[])
{
	time_t start,end;
	time(&start);

	// The following turns off sync between C and C++ streams.
	// Apparently it's faster to turn sync off, and as I don't use C streams, it's okay to turn off.
	ios_base::sync_with_stdio(false);

	parameters params(argc, argv);
	params.print_help();
	params.read_parameters();

	LOG.open((params.output_prefix + ".log").c_str());

	printLOG("\nVCFtools - " + VCFTOOLS_VERSION + "\n");
	printLOG("(C) Adam Auton 2009\n\n");

	params.print_params();

	vcf_file vcf(params.vcf_filename, params.vcf_compressed, params.chr_to_keep, params.chr_to_exclude, params.force_write_index);

	// Apply various filters as required.
	vcf.apply_filters(params);

	unsigned int N_indv = vcf.N_kept_individuals();
	unsigned int N_sites = vcf.N_kept_sites();
	printLOG("After filtering, kept " + int2str(N_indv) + " out of " + int2str(vcf.N_indv) + " Individuals\n");
	printLOG("After filtering, kept " + int2str(N_sites) + " out of a possible " + int2str(vcf.N_entries) + " Sites\n");
	if (N_sites == 0)
		error("No data left for analysis!");

	if (params.diff_file != "")
	{	// Merge files - cannot be run with other output options.
		vcf_file vcf_diff(params.diff_file, params.diff_file_compressed, params.chr_to_keep, params.chr_to_exclude, params.force_write_index);
		vcf_diff.apply_filters(params);	// Apply various filters as required.
		vcf.output_indv_in_files(params.output_prefix, vcf_diff);
		vcf.output_sites_in_files(params.output_prefix, vcf_diff);

		if (params.diff_site_discordance == true) vcf.output_discordance_by_site(params.output_prefix, vcf_diff);
		if (params.diff_discordance_matrix == true) vcf.output_discordance_matrix(params.output_prefix, vcf_diff);
		if (params.diff_indv_discordance == true) vcf.output_discordance_by_indv(params.output_prefix, vcf_diff);
		if (params.diff_switch_error == true) vcf.output_switch_error(params.output_prefix, vcf_diff);
	}

	vcf.output_INFO_for_each_site(params.output_prefix, params.INFO_to_extract);
	vcf.output_FORMAT_information(params.output_prefix, params.FORMAT_id_to_extract);
	if (params.output_indv_depth == true) vcf.output_individuals_by_mean_depth(params.output_prefix);
	if (params.output_geno_depth == true) vcf.output_genotype_depth(params.output_prefix);
	if (params.output_site_depth == true) vcf.output_site_depth(params.output_prefix, false);
	if (params.output_site_mean_depth == true) vcf.output_site_depth(params.output_prefix, true);
	if (params.output_freq == true) vcf.output_frequency(params.output_prefix, false, params.suppress_allele_output);
	if (params.output_counts == true) vcf.output_frequency(params.output_prefix, true, params.suppress_allele_output);
	if (params.plink_output == true) vcf.output_as_plink(params.output_prefix);
	if (params.plink_tped_output == true) vcf.output_as_plink_tped(params.output_prefix);
	if (params.output_HWE == true) vcf.output_hwe(params.output_prefix);
	if (params.output_SNP_density_bin_size > 0) vcf.output_SNP_density(params.output_prefix, params.output_SNP_density_bin_size);
	if (params.output_missingness == true) vcf.output_missingness(params.output_prefix);
	if (params.output_geno_rsq == true) vcf.output_genotype_r2(params.output_prefix, params.ld_snp_window_size, params.ld_bp_window_size, params.min_r2);
	if (params.output_interchromosomal_rsq == true) vcf.output_interchromosomal_genotype_r2(params.output_prefix, params.min_r2);
	if (params.output_hap_rsq == true) vcf.output_haplotype_r2(params.output_prefix, params.ld_snp_window_size, params.ld_bp_window_size, params.min_r2);
	if (params.output_het == true) vcf.output_het(params.output_prefix);
	if (params.output_site_quality == true) vcf.output_site_quality(params.output_prefix);
	if (params.output_012_matrix == true) vcf.output_as_012_matrix(params.output_prefix);
	if (params.output_as_IMPUTE == true) vcf.output_as_IMPUTE(params.output_prefix);
	if (params.output_BEAGLE_genotype_likelihoods == true) vcf.output_BEAGLE_genotype_likelihoods(params.output_prefix);
	if (params.output_as_ldhat_unphased == true) vcf.output_as_LDhat_unphased(params.output_prefix, params.chr_to_keep);
	if (params.output_as_ldhat_phased == true) vcf.output_as_LDhat_phased(params.output_prefix, params.chr_to_keep);
	if (params.output_singletons == true) vcf.output_singletons(params.output_prefix);
	if (params.output_site_pi == true) vcf.output_per_site_nucleotide_diversity(params.output_prefix);
	if (params.pi_window_size > 0) vcf.output_windowed_nucleotide_diversity(params.output_prefix, params.pi_window_size);
	if (params.output_Tajima_D_bin_size > 0) vcf.output_Tajima_D(params.output_prefix, params.output_Tajima_D_bin_size);
	if (params.output_TsTv_bin_size > 0) vcf.output_TsTv(params.output_prefix, params.output_TsTv_bin_size);
	if (params.output_TsTv_by_count) vcf.output_TsTv_by_count(params.output_prefix);
	if (params.output_TsTv_by_qual) vcf.output_TsTv_by_quality(params.output_prefix);
	if (params.recode == true) vcf.print(params.output_prefix, params.recode_INFO_to_keep, params.recode_all_INFO);
	if (params.output_filter_summary == true) vcf.output_FILTER_summary(params.output_prefix);
	if (params.output_filtered_sites == true) vcf.output_kept_and_removed_sites(params.output_prefix);
	if (params.output_LROH == true) vcf.output_LROH(params.output_prefix);
	if (params.output_relatedness == true) vcf.output_indv_relatedness(params.output_prefix);
	if (params.output_PCA == true) vcf.output_PCA(params.output_prefix, !params.PCA_no_normalisation, params.output_N_PCA_SNP_loadings);
	if (params.fst_populations.size() > 0) vcf.output_fst_version_2(params.output_prefix, params.fst_populations);

	if (params.fst_file != "")
	{
		vcf_file vcf_fst(params.fst_file, params.fst_file_compressed, params.chr_to_keep);
		vcf_fst.apply_filters(params);	// Apply various filters as required.
		vcf.output_fst(params.output_prefix, vcf_fst);
	}

	time(&end);
	double running_time = difftime(end,start);
	printLOG("Run Time = " + dbl2str_fixed(running_time, 2) + " seconds\n");
	LOG.close();
	return 0;
}
