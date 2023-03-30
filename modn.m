%trasnfer
in1='atgaacactc aaatcctggt attcgctctg attgcgatca ttccaacaaa tgcagacaaa atctgcctcg gacatcatgc cgtgtcaaac ggaaccaaag taaacacatt aactgaaaga ggagtggaag tcgtcaatgc aactgaaaca gtggaacgaa caaacacccc caggatctgc tcaaaaggga aaaggacagt tgacctcggt caatgtggac tcctggggac aatcactgga ccacctcaat gtgaccaatt cctagaattt tcggccgatt taattattga gaggcgagaa ggaagtgatg tctgttatcc tggaaaattc gtgaatgaag aagctctgag gcaaattctc agagaatcag gcggaattga caaggaaccc atgggattca catacaatgg aataagaact aatggggtga ccagtgcatg taggagatca ggatcttcat tctatgcaga aatgaaatgg ctcctgtcaa acacagataa tgctgcattc ccgcagatga ctaagtcata taaaaataca agaaaaagcc cagctataat agtatggggg atccatcatt ccgtttcaac tgcagagcaa accaagctat atgggagtgg aaacaaactg gtgacagttg ggagttctaa ttatcaacaatctttcgtac cgagtccagg agcaagacca caagttaatg gtcaatctgg aagaattgac tttcattggc taatactaaa tcccaatgat acagtcactt tcagtttcaa tggggctttc atagctccag accgtgcaag cttcctgaga ggaaaatcta tgggaatcca gagtggggta caggttgatg ccaattgtga aggggactgc tatcatagtg gagggacaat aataagtaac ttgccatttc agaacataga tagcagggca gttggaaaat gtccgagata tgttaagcaa aggagtcttc tgctggcaac agggatgaag aatgttcctg aggttccaaa gagaaaacgg actgcgagag gcctatttgg tgctatagcg ggtttcattg aaaatggatg ggaaggccta attgatggtt ggtatggttt cagacaccag aatgcacagg gagagggaac tgctgcagattacaaaagca ctcaatcggc aattgatcaa ataacaggga aattaaaccg gcttatagca aaaaccaacc aacaatttaa gttgatagac aatgaattca atgaggtaga gaagcaaatc ggtaatgtga taaattggac cagagattct ataacagaag tatggtcata caatgctgaa ctcttggtgg caatggagaa ccagcataca attgatctgg ctgattcaga aatggacaaa ctgtacgaaa gagtgaaaag acagctgaga gagaatgctg aagaagatgg cacgggttgc tttgaaatat ttcacaagtg tgatgatgac tgtatggcca gtattagaaa taacacctat gatcacagaa aatacagaga agaggcaatg caaaatagaa tacagattga cccagtcaaa ctaagcagcg gctacaaaga tgtgatactt tggtttagct tcggggcatc atgtttcata cttctagcca ttgtaatggg ccttgtcttc atatgtgtga agaatggaaa catgcggtgc actatttgta tataa';
in1=strrep(in1, ' ', '');in1=upper(in1);n=length(in1);x=zeros(n,1);
for i=1:n
    switch in1(i)
        case 'A'
            x(i)=1;
        case 'C'
            x(i)=2;
        case 'G'
            x(i)=3;
        case 'T'
            x(i)=4;
        otherwise
            error('do not input other character');
    end
end
in2='caaaaacttc ctggaaatga caatagcacg gcaacgctgt gccttgggca ccatgcagta ccaaacggaa cgatagtgaa aacaatcacg aatgaccgaa ttgaagttac taatgccact gagctggttc agaattcctc aataggtgaa atatgcgaca gtcctcatca gatccttgat ggagaaaact gcacactaat agatgctcta ttgggagacc ctcagtgtga tggctttcaa aataataaat gggacctttt tgttgaacga agcaaagcct acagcaactg ttacccttat gatgtgccgg attatgcctc ccttaggtca ctagttgcct catccggcac actggagttt aacaatgaaa gcttcaattg gactggagtc actcaaaacg gaacaagttc tgcttgcata aggagatcta atagtagttt ctttagtaga ttaaattggt tgacccactt aaacttcaaa tacccagcat tgaacgtgac tatgccaaac aatgaacaat ttgacaaatt gtacatttgg ggggttcacc acccgggtac ggacaaggac caaatcttcc tgtatgctca atcatcagga agaatcacag tatctaccaa aagaagccaa caagctgtaa tcccaaatat cggatctaga cccagaataa ggaatatccc tagcagaata agcatctatt ggacaatagt aaaaccggga gacatacttt tgattaacag cacagggaat ctaattgctc ctaggggtta cttcaaaata cgaagtggga aaagctcaat aatgagatca gatgcaccca ttggcaaatg caagtctgaa tgcatcactc caaatggaag cattcccaat gacaaaccat tccaaaatgt aaacaggatt acatacgggg cctgtcccag atatgttaag caaagcaccc tgaaattggc aacaggaatg cgaaatgtac cagagaaaca aact';
in2=strrep(in2, ' ', '');in2=upper(in2);n=length(in2);y=zeros(n,1);
for i=1:n
    switch in2(i)
        case 'A'
            y(i)=1;
        case 'C'
            y(i)=2;
        case 'G'
            y(i)=3;
        case 'T'
            y(i)=4;
        otherwise
            error('do not input other character');
    end
end
[Dist,D,k,w,rw,tw]=dtw(x,y,1);
origin=length(in1)+length(in2);
fprintf('最短距离为%d\n',Dist/origin)
%DTW
function [Dist,D,k,w,rw,tw]=dtw(r,t,pflag)
[row,M]=size(r); if (row > M) M=row; r=r'; end;
[row,N]=size(t); if (row > N) N=row; t=t'; end;
d=sqrt((repmat(r',1,N)-repmat(t,M,1)).^2);
D=zeros(size(d));
D(1,1)=d(1,1);
for m=2:M
    D(m,1)=d(m,1)+D(m-1,1);
end
for n=2:N
    D(1,n)=d(1,n)+D(1,n-1);
end
for m=2:M
    for n=2:N
        D(m,n)=d(m,n)+min(D(m-1,n),min(D(m-1,n-1),D(m,n-1)));
    end
end
 
Dist=D(M,N);
n=N;
m=M;
k=1;
w=[M N];
while ((n+m)~=2)
    if (n-1)==0
        m=m-1;
    elseif (m-1)==0
        n=n-1;
    else 
      [values,number]=min([D(m-1,n),D(m,n-1),D(m-1,n-1)]);
      switch number
      case 1
        m=m-1;
      case 2
        n=n-1;
      case 3
        m=m-1;
        n=n-1;
      end
  end
    k=k+1;
    w=[m n; w];
end
rw=r(w(:,1));
tw=t(w(:,2));
if pflag
    figure('Name','DTW - Accumulated distance matrix and optimal path', 'NumberTitle','off');
    main1=subplot('position',[0.19 0.19 0.67 0.79]);
    image(D);
    cmap = contrast(D);
    colormap(cmap);
    hold on;
    x=w(:,1); y=w(:,2);
    ind=find(x==1); x(ind)=1+0.2;
    ind=find(x==M); x(ind)=M-0.2;
    ind=find(y==1); y(ind)=1+0.2;
    ind=find(y==N); y(ind)=N-0.2;
    plot(y,x,'-w', 'LineWidth',1);
    hold off;
    axis([1 N 1 M]);
    set(main1, 'FontSize',7, 'XTickLabel','', 'YTickLabel','');
 
    colorb1=subplot('position',[0.88 0.19 0.05 0.79]);
    nticks=8;
    ticks=floor(1:(size(cmap,1)-1)/(nticks-1):size(cmap,1));
    mx=max(max(D));
    mn=min(min(D));
    ticklabels=floor(mn:(mx-mn)/(nticks-1):mx);
    colorbar(colorb1);
    set(colorb1, 'FontSize',7, 'YTick',ticks, 'YTickLabel',ticklabels);
    set(get(colorb1,'YLabel'), 'String','Distance', 'Rotation',-90, 'FontSize',7, 'VerticalAlignment','bottom');
    
    left1=subplot('position',[0.07 0.19 0.10 0.79]);
    plot(r,M:-1:1,'-b');
    set(left1, 'YTick',mod(M,10):10:M, 'YTickLabel',10*rem(M,10):-10:0)
    axis([min(r) 1.1*max(r) 1 M]);
    set(left1, 'FontSize',7);
    set(get(left1,'YLabel'), 'String','Samples', 'FontSize',7, 'Rotation',-90, 'VerticalAlignment','cap');
    set(get(left1,'XLabel'), 'String','Amp', 'FontSize',6, 'VerticalAlignment','cap');
    
    bottom1=subplot('position',[0.19 0.07 0.67 0.10]);
    plot(t,'-r');
    axis([1 N min(t) 1.1*max(t)]);
    set(bottom1, 'FontSize',7, 'YAxisLocation','right');
    set(get(bottom1,'XLabel'), 'String','Samples', 'FontSize',7, 'VerticalAlignment','middle');
    set(get(bottom1,'YLabel'), 'String','Amp', 'Rotation',-90, 'FontSize',6, 'VerticalAlignment','bottom');
    figure('Name','DTW - warped signals', 'NumberTitle','off');
    subplot(1,2,1);
    set(gca, 'FontSize',7);
    hold on;
    plot(r,'-bx');
    plot(t,':r.');
    hold off;
    axis([1 max(M,N) min(min(r),min(t)) 1.1*max(max(r),max(t))]);
    grid;
    legend('signal 1','signal 2');
    title('Original signals');
    xlabel('Samples');
    ylabel('Amplitude');
    
    subplot(1,2,2);
    set(gca, 'FontSize',7);
    hold on;
    plot(rw,'-bx');
    plot(tw,':r.');
    hold off;
    axis([1 k min(min([rw; tw])) 1.1*max(max([rw; tw]))]);
    grid;
    legend('signal 1','signal 2');
    title('Warped signals');
    xlabel('Samples');
    ylabel('Amplitude');
end
end