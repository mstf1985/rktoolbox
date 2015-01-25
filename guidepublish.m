close all
clear all
format shorte
format compact

height = 2/3;
set(0,...
'defaultfigureposition',1.5*[180 100 800 800*height],...
'defaultaxeslinewidth',1,...
'defaultaxesfontsize',16,...
'defaultlinelinewidth',2,...
'defaultpatchlinewidth',1,...
'defaultlinemarkersize',16,...
'defaulttextinterpreter','tex');

%%
publish('guide.m','imageFormat','png','stylesheet','html/mxdom2mathjax.xsl',...
    'maxWidth',600,'maxHeight',550);
movefile('html/guide.html','html/index.html')

% Read file in a string:
fid = fopen('html/index.html');
str = fscanf(fid, '%c');
fclose(fid);

% Search-replace in the string:
str = strrep(str, 'Guettel', 'G&uuml;ttel');
str = strrep(str, 'RKTOOLBOXWEBPAGE', '<p><a href="http://guettel.com/rktoolbox/rktoolbox.zip">http://guettel.com/rktoolbox/rktoolbox.zip</a></p>');

% Write back into the file;
fid = fopen('html/index.html', 'w');
fprintf(fid, '%s', str);
fclose(fid);


publish('guide.m','format','latex','outputDir','latex','stylesheet','latex/mxdom2latex.xsl');

% Read file in a string:
fid = fopen('latex/guide.tex');
str = fscanf(fid, '%c');
fclose(fid);

% Search-replace in the string:
str = strrep(str, 'Guettel', 'G{\"u}ttel');
str = strrep(str, 'RKTOOLBOXWEBPAGE', '\vspace{-1em}\begin{center}\url{http://guettel.com/rktoolbox/rktoolbox.zip}\end{center}');

% Write back into the file;
fid = fopen('latex/guide.tex', 'w');
fprintf(fid, '%s', str);
fclose(fid);


%opt.format = 'latex'; opt.outputDir = 'latex'; opt.stylesheet = 'latex/mxdom2latex.xsl'; publish('guide.m', opt)

close all
disp('Done. Guide published to html and latex.')
