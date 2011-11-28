hardcopyLattice = function(plt=NULL,
                           fileName='',
                           filePath = file.path(getwd(),'figures'),
                           eps2pdf=TRUE,
                           family='Helvetica',
                           pointsize.text = 16,
                           runQuiet=FALSE,
                           columnsPerRow = 2, #4 is probably the max==nx
                           pdf.height = 8,
                           pdf.width = 8
                          )
{
#---------------------------------------------------------------------------
#
#   This is the Multiple Panel version of hardcopyLattice.r which works
#   with the print.trellis "split" and "more" options.  All the code is the
#   same as in hardcopyLattice.r except for the multiple panels. It will
#   work with just one or with any number of trellis objects passed.
#   If one single object of class trellis is passed in 'plt' it will be
#   wrapped into a list to accommodate the multiple list code.
#
#   This routine will print a hardcopy of a lattice object. It can fail if
#   not all of the data used, say in a panel function, for adding to the
#   plot outside of the normal lattice arguments, is not available. This
#   can happen, e.g., if the lattice object is passed back from a function
#   that uses an external panel function to generate the object. In this
#   case, the extra data in that function environment are lost on the
#   function return. A much more sound way to approach returning lattice
#   objects, and one that will insure that hardcopyLattice does not fail,
#   is to always define the panel function with the environment (i.e.,
#   function) creating the plot. Then R's closure mechanism will return a
#   copy of the function environment along with the lattice object,
#   containing any extra data defined within the function and possibly
#   used in the creation of the plot. This is an important R programming
#   construct/nuance that should be heeded. 
#
#   Arguments...
#     plt = a list of lattice plot object of class 'trellis' OR a single
#           object of class trellis that will be wrapped into a list for
#           processing
#     fileName = '' for no hardcopy, a name with png, eps, or pdf extension
#                otherwise.
#     filePath = the complete path to fileName
#     eps2pdf = T: convert eps to pdf; F output as eps
#     family = the font family
#     pointsize.text = point size for text (see lattice comment below)
#     runQuiet = T: no info cat-ed; F: tell what was done
#     columnsPerRow = number of columns for each row in output matrix
#     pdf.height = in inches; note that this can affect the way panels are
#                  presented. For example, if numPlts=4 and pdf.height=8,
#                  with 2 time periods (2 panels), they are side-by-side;
#                  but if pdf.height=10, they are one atop the other in each
#                  of the 4 split displays.
#     pdf.width = width of a pdf output; making this larger that pdf.height
#                 is one way to get "landscape" in pdf output
#
#   Returns...
#     Nothing, invisibly...is that a kind of double negative?
#
#Author...									Date: 9-Jan-2008
#	Jeffrey H. Gove
#	USDA Forest Service
#	Northern Research Station
#	271 Mast Road
#	Durham, NH 03824
#	jhgove@unh.edu
#	phone: 603-868-7667	fax: 603-868-7604
#---------------------------------------------------------------------------
#
#   a couple quick checks...
#
    if(class(plt)=='trellis')  #make a list of plots passed singly as trellis objects
      plt = list(plt)
    if(!all(sapply(plt, class)=='trellis'))
      stop('Individual list objects passed not all of class trellis!')
    if(nchar(fileName) == 0)
      stop('No output file name requested!')

    numPlts = length(plt)
    if(numPlts == 0)
      stop('No plots to work with!')

#
#   lattice does not allow specifying the 'pointsize' option in the output
#   devices directly, we must change it here...
#
    if(!is.na(pointsize.text))
      for(i in 1:numPlts)
        plt[[i]] = update(plt[[i]], par.settings=list(fontsize=list(text=pointsize.text)))


    columnsPerRow = ifelse(numPlts == 1, 1, columnsPerRow) #just plotMarginalPDFs.r()
    n.rows = ceiling(numPlts/columnsPerRow) #these are vertical==ny
    idxPlts = 0   

#
#   little function to do the actual printing over all trellis objects...
#
    printMe = function() {
                 for(i in 1:n.rows) {
                   for(j in 1:columnsPerRow) {
                     idxPlts = idxPlts + 1
                     if(idxPlts < numPlts) #i.e., n.states
                       more = TRUE
                     else
                       more = FALSE
                     print(plt[[idxPlts]], split=c(j,i, columnsPerRow, n.rows), more=more)
                     if(idxPlts == numPlts)
                       break
                   } #columns
                 }  #rows             
        }#printME

    
#
#   png, tiff, eps & pdf are supported...
#
    devCur = dev.cur()
    if(!runQuiet)
      cat('\nPlotting to file...\n')
    fn = unlist(strsplit(fileName, "\\."))
    switch(fn[2],
        png = {
             fname.png = paste(fn[1], '.png', sep='')
             fileName.png <- file.path(filePath, fname.png)
             if(!runQuiet)
               cat(fileName.png, '\n')
             trellis.device('png', file= fileName.png)
             printMe()
             dev.off()
        },
        tiff = { #needs some work to get things scaled correctly...
             fname.tiff = paste(fn[1], '.tiff', sep='')
             fileName.tiff <- file.path(filePath, fname.tiff)
             if(!runQuiet)
               cat(fileName.tiff, '\n')
             trellis.device('tiff', file= fileName.tiff, antialias='subpixel',
                            width=2560, height=1920)
             plt = update(plt, par.settings = list(plot.symbol=list(cex=1.2),
                                 superpose.symbol=list(cex=rep(2,7))) #doesn't help the '+' symbols
                         )
             print(plt)
             dev.off()
        },
        eps = {
             fname.eps = paste(fn[1], '.eps', sep='')
             fileName.eps <- file.path(filePath, fname.eps)
             if(eps2pdf) {
               fname.pdf = paste(fn[1], '.pdf', sep='')
               fileName.pdf <- file.path(filePath, fname.pdf)
               if(!runQuiet)
                 cat(fileName.pdf, '\n')
             }
             else
               if(!runQuiet)
                 cat(fileName.eps, '\n')
             trellis.device('postscript', file= fileName.eps, family=family)
             printMe()
             dev.off()
             if(eps2pdf) {
               system(paste("epstopdf --nocompress",  fileName.eps, '-o', fileName.pdf))
               system(paste('rm', fileName.eps))
             }
        },
        pdf = {
               fname.pdf = paste(fn[1], '.pdf', sep='')
               fileName.pdf <- file.path(filePath, fname.pdf)
               if(!runQuiet)
                 cat(fileName.pdf, '\n')
               trellis.device('pdf', file= fileName.pdf, family=family,
                              height=pdf.height, width=pdf.width,
                              version='1.4')
               printMe()
               dev.off()
              },
        stop(paste('Illegal graphics type (',fn[2],') chosen!',sep=''))
    ) #switch
    dev.set(devCur) #preserves original device

    return(invisible())
} #hardcopyLattice
