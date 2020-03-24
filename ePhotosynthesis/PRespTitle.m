%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
%   Copyright   Xin-Guang Zhu, Yu Wang, Donald R. ORT and Stephen P. LONG
%CAS-MPG Partner Institute for Computational Biology, Shanghai Institutes for Biological Sciences, CAS, Shanghai,200031 
%China Institute of Genomic Biology and Department of Plant Biology, Shanghai Institutes for Biological Sciences, CAS, Shanghai,200031 
%University of Illinois at Urbana Champaign
%Global Change and Photosynthesis Research Unit, USDA/ARS, 1406 Institute of Genomic Biology, Urbana, IL 61801, USA.
 
%   This file is part of e-photosynthesis.
 
%    e-photosynthesis is free software; you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation; 
 
%    e-photosynthesis is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
 
%    You should have received a copy of the GNU General Public License (GPL)
%    along with this program.  If not, see <http://www.gnu.org/licenses/>.
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function suc = PRespTitle(m,p,n)
      
       if n==1
            if m ==1
                title('Glycerate');
            elseif m==2
                title('Glycolate');
            elseif m==3
                title('PGA');
            elseif m==4
                title('PGCA');
            elseif m==5
                title('C-glycolate');
            elseif m==6
                title('C-glyoxylate');
            elseif m==7
                title('C-Serine');
            elseif m==8
                title('C-Glycine');
            elseif m==9
                title('C-HPR');
            else m==10
                title('C-glycerate');
           end

       elseif n ==2
           
            if m ==1
                title('Rubisco oxygenation');
            elseif m==2
                title('Pglycolate phosphatase');
            elseif m==3
                title('Glycerate kinase');
            elseif m==4
                title('Pglycolate oxidase');
            elseif m==5
                title('Serine Glyoxylate transaminase');
            elseif m==6
                title('HOpyruvate oxidase');
            elseif m==7
                title('Glu glyoxylate transaminase');
            elseif m==8
                title('Glycine decarboxylase');
            elseif m==9
                title('Glycerate uptake');
            else m==10
                title('Glycolate export');
           end
       end
       
       suc = 1;