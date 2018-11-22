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
      'subset','s',2,'character',
      'method','n',2,'character',
    'distance','d',2,'character',
     'kingdom','k',2,'character',
    'plottype','e',2,'numeric',
    'category','g',2,'numeric',
         'log','l',2,'character',
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
### select a kingdom for phyloseq plot (e.g., "phylum")
#kingdom_str<-colnames(tax_table)[options$kingdom]
kingdom_str<-options$kingdom
distance<-options$distance
plottype<-options$plottype

### prepare the directory and file name
pdffile <- gsub("[ ]+", "", paste(options$outdir,"/pdffile.pdf"))
pngfile_nmds <- gsub("[ ]+", "", paste(options$outdir,"/nmds.png"))
pngfile_nmds_facet <- gsub("[ ]+", "", paste(options$outdir,"/nmds_facet.png"))
htmlfile <- gsub("[ ]+", "", paste(options$htmlfile))
output_summary <- gsub("[ ]+", "", paste(options$log))

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
    p1<-plot_ordination(phyloseq_obj,phyloseq_ord,type="taxa",color=kingdom_str,title="taxa") + theme(legend.position="bottom",legend.box="vertical",legend.direction="horizontal")
    print(p1)
    p2<-plot_ordination(phyloseq_obj,phyloseq_ord,type="taxa",color=kingdom_str,title="taxa") + facet_wrap(formula(paste('~',kingdom_str)),3) +  theme(legend.position="bottom",legend.box="vertical",legend.direction="horizontal")
    print(p2)
    garbage<-dev.off();

    #png('nmds.png')
    bitmap(pngfile_nmds,"png16m",res=120)
    p3<-plot_ordination(phyloseq_obj,phyloseq_ord,type="taxa",color=kingdom_str,title="taxa") +  theme(legend.position="bottom",legend.box="vertical",legend.direction="horizontal")
    print(p3)
    garbage<-dev.off()

    #png('nmds_facet.png')
    bitmap(pngfile_nmds_facet,"png16m",res=120)
    p4<-plot_ordination(phyloseq_obj,phyloseq_ord,type="taxa",color=kingdom_str,title="taxa") + facet_wrap(formula(paste('~',kingdom_str)),3) +  theme(legend.position="bottom",legend.box="vertical",legend.direction="horizontal")
    print(p4)
    garbage<-dev.off()

    create_HTML_1(htmlfile)
}

create_SAMPLE_PDF<-function(pdf_file,phyloseq_obj,phyloseq_ord,htmlfile,pngfile_nmds,category_type){
    pdf(pdf_file);
    p <- plot_ordination(phyloseq_obj, phyloseq_ord, type="samples", color=category_type)
    p <- p + geom_point(aes(fill=category_type)) + geom_point(size=5) + ggtitle(paste("Samples - Stress value",formatC(phyloseq_ord$stress,digits=4,format="f"),sep=":")) +  theme(legend.position="bottom",legend.box="vertical",legend.direction="horizontal")
    print(p)
    garbage<-dev.off();

    #png('nmds.png')
    bitmap(pngfile_nmds,"png16m",res=120)
    p1 <- plot_ordination(phyloseq_obj, phyloseq_ord, type="samples", color=category_type)
    p1 <- p1 + geom_point(aes(fill=category_type)) + geom_point(size=5) + ggtitle(paste("Samples - Stress value",formatC(phyloseq_ord$stress,digits=4,format="f"),sep=":")) +  theme(legend.position="bottom",legend.box="vertical",legend.direction="horizontal")
    print(p1)
    garbage<-dev.off();

    create_HTML_2(htmlfile)
}

create_BIPLOT_PDF<-function(pdf_file,phyloseq_obj,phyloseq_ord,kingdom_str,htmlfile,pngfile_nmds,category_type){
    pdf(pdf_file);
    print(category_type)
    p_biplot <- plot_ordination(phyloseq_obj, phyloseq_ord, type="biplot", color=category_type, shape=kingdom_str,title="BIPLOT") +  theme(legend.position="bottom",legend.box="vertical",legend.direction="horizontal")
    print(p_biplot)
    garbage<-dev.off();

    bitmap(pngfile_nmds,"png16m",res=120)
    p_biplot_png <- plot_ordination(phyloseq_obj, phyloseq_ord, type="biplot", color=category_type, shape=kingdom_str,title="BIPLOT") +  theme(legend.position="bottom",legend.box="vertical",legend.direction="horizontal")
    print(p_biplot_png)
    garbage<-dev.off();

    create_HTML_2(htmlfile)
}

create_SPLITPLOT_PDF<-function(pdf_file,phyloseq_obj,phyloseq_ord,kingdom_str,htmlfile,pngfile_nmds,category_type){
   pdf(pdf_file,width=10, height=6);
   split_plot <- plot_ordination(phyloseq_obj, phyloseq_ord, type="split", color=kingdom_str, shape=kingdom_str, label=category_type, title="SPLIT PLOT")
   split_plot <- split_plot + theme(plot.margin = unit(c(12,18,12,18),"pt"),legend.position="bottom",legend.box="vertical",legend.direction="horizontal")
   print(split_plot)
   garbage<-dev.off();

   bitmap(pngfile_nmds,"png16m", res=120)
   split_plot <- plot_ordination(phyloseq_obj, phyloseq_ord, type="split", color=kingdom_str, shape=kingdom_str, label=category_type, title="SPLIT PLOT") 
   split_plot <- split_plot + theme(plot.margin = unit(c(12,18,12,18),"pt"),legend.position="bottom",legend.box="vertical",legend.direction="horizontal")
   print(split_plot)
   garbage<-dev.off();
   create_HTML_2(htmlfile)
}

create_HTML_1<-function(htmlfile){
    htmlfile_handle <- file(htmlfile)
    html_output = c('<html><body>',
                    '<table align="center">',
                    '<tr>',
                    '<td valign="middle" style="vertical-align:middle;">',
                        '<a href="pdffile.pdf"><img src="nmds.png" width="800" height="800"/></a>',
                    '</td>',
                    '</tr>',
                    '<tr>',
                    '<td valign="middle" style="vertical-align:middle;">',
                    '<a href="pdffile.pdf"><img src="nmds_facet.png" width="800" height="800"/></a>',
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
                    '<table align="center">',
                    '<tr>',
                    '<td valign="middle" style="vertical-align:middle;">',
                        '<a href="pdffile.pdf"><img src="nmds.png" width="800" height="800"/></a>',
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

    category_input = get_variable(physeq, category_type) %in% category_option
    sample_data(physeq)$category_input <- factor(category_input)

#   compute distance matrix
    physeq_ord<-ordinate(physeq,method,distance)
    
#   get column sum
    sum_table<-data.frame(column_sum=as.matrix(colSums(otu_table(physeq))))

    rowname_table<-data.frame(sample=rownames(sum_table))

    output_table<-as.data.frame(cbind(rowname_table,sum_table))

    output_table<-output_table[order(output_table$column_sum),]
    
#   Reformat distance matrix
    distance_matrix<-as.data.frame(physeq_ord$points)
    distance_matrix<-cbind(sample=rownames(distance_matrix),distance_matrix)

    sink(output_summary)
    cat('--------------------------------------')
    cat('\n')
    cat('Stress value')
    cat('\n')
    cat(formatC(physeq_ord$stress,digits=4,format="f"))
    cat('\n')
    cat('--------------------------------------')
    cat('\n')
    cat('Sample - Column Sum')
    cat('\n')
    cat('--------------------------------------')
    cat('\n')
    write.table(output_table,row.names=F,quote=F)
    cat('\n')
    cat('--------------------------------------')
    cat('\n')
    cat('Distance Matrix')
    cat('\n')
    cat('--------------------------------------')
    cat('\n')
    write.table(distance_matrix,row.names=F,quote=F)
    cat('\n')
    cat('--------------------------------------')
    sink()

if(plottype == 1){
    #kingdom_str = colnames(tax_table)[2]
    create_OTU_PDF(pdffile,physeq,physeq_ord,kingdom_str,htmlfile,pngfile_nmds,pngfile_nmds_facet) 
}else if(plottype == 2){
    create_SAMPLE_PDF(pdffile,physeq,physeq_ord,htmlfile,pngfile_nmds,category_type)
}else if(plottype == 3){
    create_BIPLOT_PDF(pdffile,physeq,physeq_ord,kingdom_str,htmlfile,pngfile_nmds,category_type)
}else{
    create_SPLITPLOT_PDF(pdffile,physeq,physeq_ord,kingdom_str,htmlfile,pngfile_nmds,category_type)
}

