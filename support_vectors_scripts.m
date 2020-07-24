svm = fitcsvm( ch1_ch2_prebleach_F, (fret_by_particle>0.7) )
%%
mydata = ch1_ch2_prebleach_F;

sv = svm.SupportVectors
figure; plot( sv(:,1), sv(:,2), 'o' )

output = polyfit( sv(:,1), sv(:,2) , 1 );
[ xmin,xmax ] = deal( -output(2)/output(1)*100, max(mydata(:,1)) )
vertices = [ xmin, polyval(output,xmin) ; xmax, polyval(output,xmin); xmax, polyval(output,xmax) ]

in_poly = inpolygon( mydata(:,1), mydata(:,2), vertices(:,1), vertices(:,2) );

figure; plot( sv(:,1), sv(:,2), 'o' ); lsline
hold on;  plot( ch1_ch2_prebleach_F( in_poly, 1 ), ch1_ch2_prebleach_F( in_poly, 2 ), 'o' )
hold on;  plot( ch1_ch2_prebleach_F( ~in_poly, 1 ), ch1_ch2_prebleach_F( ~in_poly, 2 ), 'ko' )

[x1,x2] = deal( ch1_ch2_prebleach_F( ~in_poly, 1 ), ch1_ch2_prebleach_F( ~in_poly, 2 ) );
GMModel = fitgmdist ( [x1, x2], 1 );

gmPDF = @(x1,x2)reshape(pdf(GMModel,[x1(:) x2(:)]),size(x1));
g = gca;
fcontour(gmPDF,[g.XLim g.YLim])

figure; scatter( ch1_ch2_prebleach_F(:,1), ch1_ch2_prebleach_F(:,2), 20, gmPDF( ch1_ch2_prebleach_F(:,1), ch1_ch2_prebleach_F(:,2) ) );
%plot( vertices(:,1), vertices(:,2), 'r--' )

%%

figure; plot( sv(:,1), sv(:,2), 'o' ); lsline; hold on;  plot( ch1_ch2_prebleach_F(find( svm.W == min( svm.W ) ),1),ch1_ch2_prebleach_F(find( svm.W == min( svm.W ) ),2),'.' );
threshold_ = max( diff( sort(gmPDF( ch1_ch2_prebleach_F(:,1), ch1_ch2_prebleach_F(:,2) )) ) );

gm_ = gmPDF( ch1_ch2_prebleach_F(:,1), ch1_ch2_prebleach_F(:,2) );

find(gm_>threshold_)