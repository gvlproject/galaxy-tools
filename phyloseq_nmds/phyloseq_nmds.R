library('getopt')
library('ape')
suppressPackageStartupMessages(library('phyloseq'))
library(biomformat)
library(plyr)
Sys.setenv("DISPLAY"=":1")
library("ggplot2")
suppressPackageStartupMessages(library("doParallel"))
ncores = ceiling(detectCores() * 0.8)
registerDoParallel(cores=ncores)

options(warn=-1)

theme_set(theme_bw())

#http://saml.rilspace.com/creating-a-galaxy-tool-for-r-scripts-that-output-images-and-pdfs
#http://joey711.github.io/phyloseq-demo/phyloseq-demo.html
option_specification = matrix(c(
   'otu_table','o',2,'character',
   'tax_table','t',2,'character',
  'meta_table','m',2,'character',
        'biom','b',2,'character',
      'filter','f',2,'numeric',
      'subset','s',2,'character',
      'method','n',2,'character',
    'distance','d',2,'character',
     'kingdom','k',2,'character',
    'plottype','e',2,'numeric',
      'cutoff','c',2,'numeric',
    'category','g',2,'numeric',
        'keep','p',2,'numeric',
      'outdir','r',2,'character',
    'htmlfile','h',2,'character'
),byrow=TRUE,ncol=4);


options <- getopt(option_specification);
options(bitmapType="cairo")

if (!is.null(options$outdir)) {
  # Create the directory
  dir.create(options$outdir,FALSE)
}



method<-options$method
cutoff_value<-options$cutoff
### select a kingdom for phyloseq plot (e.g., "phylum")
#kingdom_str<-colnames(tax_table)[options$kingdom]
kingdom_str<-options$kingdom
distance<-options$distance
keep<-options$keep
plottype<-options$plottype
filter<-options$filter

### prepare the directory and file name
pdffile <- gsub("[ ]+", "", paste(options$outdir,"/pdffile.pdf"))
pngfile_nmds <- gsub("[ ]+", "", paste(options$outdir,"/nmds.png"))
pngfile_nmds_facet <- gsub("[ ]+", "", paste(options$outdir,"/nmds_facet.png"))
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


create_OTU_PDF<-function(pdf_file,phyloseq_obj,phyloseq_ord,kingdom_str,htmlfile,pngfile_nmds,pngfile_nmds_facet){
    pdf(pdf_file);
    p1<-plot_ordination(physeq,physeq_ord,type="taxa",color="Phylum",title="taxa")
    print(p1)
    p2<-plot_ordination(physeq,physeq_ord,type="taxa",color="Phylum",title="taxa") + facet_wrap(formula(paste('~',kingdom_str)),3)
    print(p2)
    garbage<-dev.off();

    #png('nmds.png')
    bitmap(pngfile_nmds,"png16m")
    p3<-plot_ordination(physeq,physeq_ord,type="taxa",color="Phylum",title="taxa")
    print(p3)
    garbage<-dev.off()

    #png('nmds_facet.png')
    bitmap(pngfile_nmds_facet,"png16m")
    p4<-plot_ordination(physeq,physeq_ord,type="taxa",color="Phylum",title="taxa") + facet_wrap(formula(paste('~',kingdom_str)),3)
    print(p4)
    garbage<-dev.off()

    create_HTML_1(htmlfile)
}

create_SAMPLE_PDF<-function(pdf_file,phyloseq_obj,phyloseq_ord,htmlfile,pngfile_nmds,category_type){
    pdf(pdf_file);
    p <- plot_ordination(physeq_filter, physeq_ord, type="samples", color=category_type)
    p <- p + geom_point(aes(fill=category_type)) + geom_point(size=5) + ggtitle("samples")
    print(p)
    garbage<-dev.off();

    #png('nmds.png')
    bitmap(pngfile_nmds,"png16m")
    p1 <- plot_ordination(physeq_filter, physeq_ord, type="samples", color=category_type)
    p1 <- p1 + geom_point(aes(fill=category_type)) + geom_point(size=5) + ggtitle("samples")
    print(p1)
    garbage<-dev.off();

     create_HTML_2(htmlfile)
}

create_HTML_1<-function(htmlfile){
    htmlfile_handle <- file(htmlfile)
    html_output = c('<html><body>',
                    '<table align="center>',
                    '<tr>',
                    '<td valign="middle" style="vertical-align:middle;">',
                        '<a href="pdffile.pdf"><img src="nmds.png"/></a>',
                    '</td>',
                    '</tr>',
                    '<tr>',
                    '<td valign="middle" style="vertical-align:middle;">',
                    '<a href="pdffile.pdf"><img src="nmds_facet.png"/></a>',
                    '</td>',
                    '</tr>',
                    '</table>',
                    '</html></body>');
     writeLines(html_output, htmlfile_handle);
     close(htmlfile_handle);
}

create_HTML_2<-function(htmlfile){
    htmlfile_handle <- file(htmlfile)
    html_output = c('<html><body>',
                    '<table align="center>',
                    '<tr>',
                    '<td valign="middle" style="vertical-align:middle;">',
                        '<a href="pdffile.pdf"><img src="nmds.png"/></a>',
                    '</td>',
                    '</tr>',
                    '</table>',
                    '</html></body>');
    writeLines(html_output, htmlfile_handle);
    close(htmlfile_handle);
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
    category_type <- options$subset

    ### obtain the unique value in the selected column from sample table
    category_option<-unique(sample_data(physeq))[,options$subset]

}else{

    ### read the data into correct data type to create phyloseq object
    otu_table<-as.matrix(read.table(options$otu_table,header=T,sep="\t"))
    tax_table<-as.matrix(read.table(options$tax_table,header=T,sep="\t"))
    sample_table<-read.table(options$meta_table,header=T,sep="\t",stringsAsFactors=F)


    ### select column name from sample table for nmds plot
    category_type<-colnames(sample_table)[options$category]

    ### obtain the unique value in the selected column from sample table
    category_option<-unique(sample_table[,options$category])


    ### create a sample object for phyloseq
    sample_object<-sample_data(sample_table)

    ### create otu object for phyloseq
    OTU<-otu_table(otu_table, taxa_are_rows = TRUE)

    ### create tax object for phyloseq
    TAX<-tax_table(tax_table)

    ### create a phyloseq object
    physeq = phyloseq(OTU,TAX,sample_object)
}
    ### select a kingdom for phyloseq plot (e.g., "phylum")
    #kingdom_str<-colnames(tax_table)[options$kingdom]
    #kingdom_str<-options$kingdom
    
    
    ### filtering OTUs based on cutoff value (e.g., 5)
    #physeq_temp =genefilter_sample(physeq, filterfun_sample(function(x) x > cutoff_value), A=0.1*nsamples(physeq))
    physeq_temp =genefilter_sample(physeq, filterfun_sample(function(x) x > cutoff_value), A=filter*nsamples(physeq))
    
    ### phyloseq object after filtered
    physeq_filter = prune_taxa(physeq_temp, physeq)
    
    ## Transform to even sampling depth
    physeq_filter = transform_sample_counts(physeq_filter, function(x) 1E6 * x/sum(x))



if(plottype == 1){

    # kingdom_str <- as.numeric(kingdom_str)
    ## Keep only the most abundant five phyla

    phylum.sum = tapply(taxa_sums(physeq_filter), tax_table(physeq_filter)[, kingdom_str], sum,na.rm=TRUE)
    ### number of most abundance phyla to keep
    topphyla = names(sort(phylum.sum, TRUE))[1:keep]

    physeq_filter = prune_taxa((tax_table(physeq_filter)[, kingdom_str] %in% topphyla), physeq_filter)

}else{
 
    otu_table(physeq_filter)<-otu_table(physeq_filter)[,colSums(otu_table(physeq_filter)) > 1]

}
    ### select category to plot NMDS
    category_input = get_variable(physeq_filter, category_type) %in% category_option
    sample_data(physeq_filter)$category_input <- factor(category_input)

    physeq_ord<-ordinate(physeq_filter,method,distance)

if(plottype == 1){
    #kingdom_str = colnames(tax_table)[2]
    create_OTU_PDF(pdffile,physeq_filter,physeq_ord,kingdom_str,htmlfile,pngfile_nmds,pngfile_nmds_facet) 
}else{
    create_SAMPLE_PDF(pdffile,physeq_filter,physeq_ord,htmlfile,pngfile_nmds,category_type)
}

