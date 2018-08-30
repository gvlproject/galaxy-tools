library('getopt')
library('ape')
library('ggplot2')
suppressPackageStartupMessages(library('phyloseq'))
library(biomformat)
library(plyr)
Sys.setenv("DISPLAY"=":1")
library(biomformat)
suppressPackageStartupMessages(library(metagenomeSeq))
suppressPackageStartupMessages(library("doParallel"))
ncores = ceiling(detectCores() * 0.8)
registerDoParallel(cores=ncores)

options(warn=-1)

theme_set(theme_bw())
## ggplot 
# http://r-statistics.co/Top50-Ggplot2-Visualizations-MasterList-R-Code.html
# https://gist.github.com/Mikeyj/5429538
# http://microbiome-tutorials.readthedocs.io/en/latest/_static/Composition.html
# https://rstudio-pubs-static.s3.amazonaws.com/268156_d3ea37937f4f4469839ab6fa2c483842.html#otus_that_differ_by (stacked bar plot)

## color
## http://www.cookbook-r.com/Graphs/Colors_(ggplot2)/ 

#http://saml.rilspace.com/creating-a-galaxy-tool-for-r-scripts-that-output-images-and-pdfs
#http://joey711.github.io/phyloseq-demo/phyloseq-demo.html
option_specification = matrix(c(
   'otu_table','o',2,'character',
   'tax_table','t',2,'character',
  'meta_table','m',2,'character',
        'biom','b',2,'character',
      'filter','f',2,'numeric',
     'kingdom','k',2,'character',
      'cutoff','c',2,'numeric',
        'keep','p',2,'numeric',
     'outbiom','h',2,'character',
      'outdir','d',2,'character',
    'htmlfile','w',2,'character'
),byrow=TRUE,ncol=4);


options <- getopt(option_specification);
options(bitmapType="cairo")

if (!is.null(options$outdir)) {
  # Create the directory
  dir.create(options$outdir,FALSE)
}



cutoff_value<-options$cutoff
### select a kingdom for phyloseq plot (e.g., "phylum")
#kingdom_str<-colnames(tax_table)[options$kingdom]
kingdom_str<-options$kingdom
keep<-options$keep
filter<-options$filter

### prepare the directory and file name
pdffile <- gsub("[ ]+", "", paste(options$outdir,"/pdffile.pdf"))
pngfile_before_filtering <- gsub("[ ]+", "", paste(options$outdir,"/barplot_before_filtering.png"))
pngfile_after_filtering <- gsub("[ ]+", "", paste(options$outdir,"/barplot_after_filtering.png"))
pngfile_pre_phyla_filtering <- gsub("[ ]+", "", paste(options$outdir,"/barplot_before_phyla_filtering.png"))
pngfile_post_phyla_filtering<- gsub("[ ]+", "", paste(options$outdir,"/barplot_after_phyla_filtering.png"))
htmlfile <- gsub("[ ]+", "", paste(options$htmlfile))

### This function accepts different two different type of BIOM file format
readBIOM<-function(inBiom){
	tryCatch({
		phyloseq_obj<-import_biom(inBiom,parallel=TRUE)
		return(phyloseq_obj)
		},
		error=function(e){
		biom_obj<-read_biom(inBiom)

		otu_matrix = as(biom_data(biom_obj), "matrix")
		OTU_TABLE = otu_table(otu_matrix, taxa_are_rows=TRUE)
		
		taxonomy_matrix = as.matrix(observation_metadata(biom_obj), rownames.force=TRUE)	
		TAXONOMY_TABLE = tax_table(taxonomy_matrix)	
		
		metadata.temp<-sample_metadata(biom_obj)
		METADATA_TABLE<-plyr::ldply(metadata.temp, rbind)
        rownames(METADATA_TABLE)<-as.character(METADATA_TABLE$.id)

		phyloseq_obj = phyloseq(OTU_TABLE, TAXONOMY_TABLE,sample_data(METADATA_TABLE))
		return(phyloseq_obj)
		}
	)
}


create_PDF<-function(pdf_file,OTU_DATAFRAME_BEFORE_FILTERING,OTU_DATAFRAME_AFTER_FILTERING,physeq_pre_phyla_filtering,physeq_post_phyla_filtering,kingdom_str,htmlfile,pngfile_before_filtering,pngfile_after_filtering,pngfile_pre_phyla_filtering,pngfile_post_phyla_filtering){
    pdf(pdf_file);
    barplot_before_filtering<-ggplot(OTU_DATAFRAME_BEFORE_FILTERING,aes(rownames(OTU_DATAFRAME_BEFORE_FILTERING))) + 
                                  geom_bar(aes(weight=Abundance),fill="tomato3") +
                                  labs(title="Sample Depth Bar Chart",subtitle="Sample Vs Abundance (Before Filtering)",caption="source: Input Biom") +
                                  xlab("Sample") +
                                  ylab("Abundance") +
                                  theme(axis.text.x=element_text(angle=65,vjust=0.6))
    print(barplot_before_filtering)

    barplot_after_filtering<-ggplot(OTU_DATAFRAME_AFTER_FILTERING,aes(rownames(OTU_DATAFRAME_AFTER_FILTERING))) + 
                                  geom_bar(aes(weight=Abundance),fill="blue") + 
                                  labs(title="Sample Depth Bar Chart",subtitle="Sample Vs Abundance (After Filtering)", caption="source: Input Biom") +
                                  xlab("Sample") +
                                  ylab("Abundance") +
                                  theme(axis.text.x=element_text(angle=65,vjust=0.6))
    print(barplot_after_filtering)

    barplot_pre_phyla_filtering<-plot_bar(physeq_pre_phyla_filtering, x=colnames(sample_data(physeq_pre_phyla_filtering))[1], fill=kingdom_str) + 
                                 geom_bar(stat="identity", position="stack") +
                                 labs(title=paste("Sample Depth Bar Chart",kingdom_str,sep=":"),subtitle="Sample Vs Abundance (Pre phyla Filtering)",caption="source: Input Biom") +
                                 xlab("Sample") +
                                 ylab("Abundance") +
                                 theme(axis.text.x=element_text(angle=90,vjust=0.6)) +
                                 scale_fill_hue()
    print(barplot_pre_phyla_filtering)

    barplot_post_phyla_filtering<-plot_bar(physeq_post_phyla_filtering, x=colnames(sample_data(physeq_post_phyla_filtering))[1], fill=kingdom_str) +
                                  geom_bar(stat="identity", position="stack") +
                                  labs(title=paste("Sample Depth Bar Chart",kingdom_str,sep=":"),subtitle="Sample Vs Abundance (Post phyla Filtering)",caption="source: Input Biom") +
                                  xlab("Sample") +
                                  ylab("Abundance") +
                                  theme(axis.text.x=element_text(angle=90,vjust=0.6)) +
                                  scale_fill_hue()
    print(barplot_post_phyla_filtering)
    garbage<-dev.off();

    #png('barplot_before_filtering.png')
    bitmap(pngfile_before_filtering,"png16m")
    barplot_before_filtering_png<-ggplot(OTU_DATAFRAME_BEFORE_FILTERING,aes(rownames(OTU_DATAFRAME_BEFORE_FILTERING))) +
                                  geom_bar(aes(weight=Abundance),fill="tomato3") +
                                  labs(title="Sample Depth Bar Chart",subtitle="Sample Vs Abundance (Before Filtering)",caption="source: Input Biom") +
                                  xlab("Sample") +
                                  ylab("Abundance") +
                                  theme(axis.text.x=element_text(angle=65,vjust=0.6))
    print(barplot_before_filtering_png)
    garbage<-dev.off()

    #png('barplot_after_filtering.png')
    bitmap(pngfile_after_filtering,"png16m")
    barplot_after_filtering_png<-ggplot(OTU_DATAFRAME_AFTER_FILTERING,aes(rownames(OTU_DATAFRAME_AFTER_FILTERING))) +
                                 geom_bar(aes(weight=Abundance),fill="blue") +
                                 labs(title="Sample Depth Bar Chart",subtitle="Sample Vs Abundance (After Filtering)", caption="source: Input Biom") +
                                 xlab("Sample") +
                                 ylab("Abundance") +
                                 theme(axis.text.x=element_text(angle=65,vjust=0.6))
    print(barplot_after_filtering_png)
    garbage<-dev.off()

    #png('barplot_pre_phyla_filtering.png')
    bitmap(pngfile_pre_phyla_filtering,"png16m")
    print(sample_data(physeq_pre_phyla_filtering))
    barplot_pre_phyla_filtering<-plot_bar(physeq_pre_phyla_filtering, x=colnames(sample_data(physeq_pre_phyla_filtering))[1], fill=kingdom_str) +
                                 geom_bar(stat="identity", position="stack") +
                                 labs(title=paste("Sample Depth Bar Chart",kingdom_str,sep=":"),subtitle="Sample Vs Abundance (Pre Phyla Filtering)",caption="source: Input Biom") +
                                 xlab("Sample") +
                                 ylab("Abundance") +
                                 theme(axis.text.x=element_text(angle=90,vjust=0.6)) +
                                 scale_fill_hue()
    print(barplot_pre_phyla_filtering)
    garbage<-dev.off()

    #png('barplot_post_phyla_filtering.png')
    bitmap(pngfile_post_phyla_filtering,"png16m")
    barplot_post_phyla_filtering<-plot_bar(physeq_post_phyla_filtering, x=colnames(sample_data(physeq_post_phyla_filtering))[1], fill=kingdom_str) +
                                  geom_bar(stat="identity", position="stack") +
                                  labs(title=paste("Sample Depth Bar Chart",kingdom_str,sep=":"),subtitle="Sample Vs Abundance (Post Phyla Filtering)",caption="source: Input Biom") +
                                  xlab("Sample") +
                                  ylab("Abundance") +
                                  theme(axis.text.x=element_text(angle=90,vjust=0.6)) +
                                  scale_fill_hue()
    print(barplot_post_phyla_filtering)
    garbage<-dev.off()

    create_HTML(htmlfile)
}

create_HTML<-function(htmlfile){
    htmlfile_handle <- file(htmlfile)
    html_output = c('<html><body>',
                    '<table align="center>',
                    '<tr>',
                    '<td valign="middle" style="vertical-align:middle;">',
                        '<a href="pdffile.pdf"><img src="barplot_before_filtering.png"/></a>',
                    '</td>',
                    '</tr>',
                    '<tr>',
                    '<td valign="middle" style="vertical-align:middle;">',
                    '<a href="pdffile.pdf"><img src="barplot_after_filtering.png"/></a>',
                    '</td>',
                    '</tr>',
                    '<tr>',
                    '<td valign="middle" style="vertical-align:middle;">',
                    '<a href="pdffile.pdf"><img src="barplot_before_phyla_filtering.png"/></a>',
                    '</td>',
                    '</tr>',
                    '<tr>',
                    '<td valign="middle" style="vertical-align:middle;">',
                    '<a href="pdffile.pdf"><img src="barplot_after_phyla_filtering.png"/></a>',
                    '</td>',
                    '</tr>',
                    '</table>',
                    '</html></body>');
     writeLines(html_output, htmlfile_handle);
     close(htmlfile_handle);
}

convert_phyloseq_otutable_to_dataframe<-function(physeq_obj){
    temp.df<-data.frame(otu_table(physeq_obj))   
    temp.df.counts<-as.data.frame(colSums(temp.df))
    colnames(temp.df.counts)<-"Abundance"
    print(temp.df.counts)
    return(temp.df.counts)
}

if(!is.null(options$biom)){
    
    #physeq<-import_biom(options$biom)
    physeq<-readBIOM(options$biom)
    
    if(length(rank_names(physeq)) == 8){
        tax_table(physeq) <- tax_table(physeq)[,-1]
        colnames(tax_table(physeq)) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
    } else {
        colnames(tax_table(physeq)) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
    }
    
    ### select column name from sample table for nmds plot
    ## which(colnames(sample_data(biom)) == "vegetation_type_id")
    #category_type<-colnames(sample_data(physeq))[options$subset]
    #category_type <- options$subset

    ### obtain the unique value in the selected column from sample table
    #category_option<-unique(sample_data(physeq))[,options$subset]

}else{

    ### read the data into correct data type to create phyloseq object
    otu_table<-as.matrix(read.table(options$otu_table,header=T,sep="\t"))
    tax_table<-as.matrix(read.table(options$tax_table,header=T,sep="\t"))
    sample_table<-read.table(options$meta_table,header=T,sep="\t",stringsAsFactors=F)


    ### select column name from sample table for nmds plot
    #category_type<-colnames(sample_table)[options$category]

    ### obtain the unique value in the selected column from sample table
    #category_option<-unique(sample_table[,options$category])


    ### create a sample object for phyloseq
    sample_object<-sample_data(sample_table)

    ### create otu object for phyloseq
    OTU<-otu_table(otu_table, taxa_are_rows = TRUE)

    ### create tax object for phyloseq
    TAX<-tax_table(tax_table)

    ### create a phyloseq object
    physeq = phyloseq(OTU,TAX,sample_object)
}
    ### make the first column to be the sample ID in the phyloseq object
   
    firstColumn = sample_data(physeq)[,1]
    row_names = rownames(sample_data(physeq))
    check = all( firstColumn == row_names)
    if(!check){
        sample_data(physeq) <- cbind(SampleID= rownames(sample_data(physeq)),sample_data(physeq))
    }

    
    ### extract otu table from phyloseq object
    before_filtering_dataframe_sampleCounts<-convert_phyloseq_otutable_to_dataframe(physeq)
    
    ### filtering OTUs based on cutoff value (e.g., 5)
    #physeq_temp =genefilter_sample(physeq, filterfun_sample(function(x) x > cutoff_value), A=0.1*nsamples(physeq))
    physeq_temp =genefilter_sample(physeq, filterfun_sample(function(x) x > cutoff_value), A=filter*nsamples(physeq))
    
    ### phyloseq object after filtered
    physeq_filter = prune_taxa(physeq_temp, physeq)
    
  
    
    ## Transform to even sampling depth
    #physeq_filter = transform_sample_counts(physeq_filter, function(x) 1E6 * x/sum(x))

    #after_filtering.dataframe<-data.frame(otu_table(physeq_filter))
    #after_filtering_dataframe_sampleCounts<-as.data.frame(colSums(after_filtering.dataframe))
    #colnames(after_filtering_dataframe_sampleCounts)<-"Abundance"
    after_filtering_dataframe_sampleCounts<-convert_phyloseq_otutable_to_dataframe(physeq_filter)
   
#    create_PDF(pdffile,before_filtering_dataframe_sampleCounts,after_filtering_dataframe_sampleCounts,htmlfile,pngfile_before_filtering,pngfile_after_filtering)

    # kingdom_str <- as.numeric(kingdom_str)
    ## Keep only the most abundant five phyla

    ### Phyla - Pre transformation (Transform to even sampling depth)
    
    #physeq_filter_pre_transform = tapply(taxa_sums(physeq_filter), tax_table(physeq_filter)[, kingdom_str], sum,na.rm=TRUE)
    
    phylum.sum_pre_transform= tapply(taxa_sums(physeq_filter), tax_table(physeq_filter)[, kingdom_str], sum,na.rm=TRUE)
    
    topphyla_pre_transform = names(sort(phylum.sum_pre_transform, TRUE))[1:keep]
    
    physeq_filter_pre_transform = prune_taxa((tax_table(physeq_filter)[, kingdom_str] %in% topphyla_pre_transform), physeq_filter)
   
    ### Phyla - Post Transformation (Transform to even sampling depth)
    physeq_filter_post_transform = transform_sample_counts(physeq_filter, function(x) 1E6 * x/sum(x))

    phylum.sum_post_transform = tapply(taxa_sums(physeq_filter_post_transform), tax_table(physeq_filter_post_transform)[, kingdom_str], sum,na.rm=TRUE)
    ### number of most abundance phyla to keep
    topphyla_post_transform = names(sort(phylum.sum_post_transform, TRUE))[1:keep]

    physeq_filter_post_transform = prune_taxa((tax_table(physeq_filter_post_transform)[, kingdom_str] %in% topphyla_post_transform), physeq_filter_post_transform)


    create_PDF(pdffile,before_filtering_dataframe_sampleCounts,after_filtering_dataframe_sampleCounts,physeq_filter_pre_transform,physeq_filter_post_transform,kingdom_str,htmlfile,pngfile_before_filtering,pngfile_after_filtering,pngfile_pre_phyla_filtering,pngfile_post_phyla_filtering)

    ### convert phyloseq object to metagenomeSeq - preparing for BIOM output
    metagenomeSeq_obj <- phyloseq_to_metagenomeSeq(physeq_filter_post_transform)
    metagenomeSeq_biom <- MRexperiment2biom(metagenomeSeq_obj)

    ## write biom file
    write_biom(metagenomeSeq_biom, biom_file=options$outbiom)
