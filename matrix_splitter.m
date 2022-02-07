projectdir = pwd;
dinfo = dir(fullfile(projectdir, '*.mat'));
dinfo([dinfo.isdir]) = [];     %get rid of all directories including . and ..
nfiles = length(dinfo);
for j = 1 : nfiles
	    filename = fullfile(projectdir, dinfo(j).name);
	        M = load(filename)
		    len = length(M.Problem.A)
		        newStr = split(filename,'.')
			    diags = diag(M.Problem.A);
			        new_n = strcat(newStr(1),'-diag.txt')
				    writematrix(diags',sprintf('%s',new_n{:}) ,'Delimiter','tab');
				        s = tril(M.Problem.A,-1);
					    [row, col, v] = find(s);
					        rowwise = sum(s ~= 0, 2);
						    rowptr=zeros(len+1,1);
						        [row, col, v] = find(rowwise);
							    % accumulate rowptr vector
							        for i = 1:length(row) 
									        N=zeros(len+1,1);
										        ind = row(i)+1;
											        % calculate accumulation vector
												        for j = ind:len+1
														            N(j, 1) = v(i);
															            end
																            rowptr = rowptr+N;
																	        end
																		    new_n=strcat(newStr(1),'-rowptr.txt')
																		        writematrix(rowptr,sprintf('%s',new_n{:}),'Delimiter','tab');
																			    new_n=strcat(newStr(1),'-rowwise_nonzero_count.txt')
																			        writematrix(rowwise,sprintf('%s',new_n{:}),'Delimiter','tab');
																				    new_n=strcat(newStr(1),'-row_except_diag.txt') 
																				        writematrix(row,sprintf('%s',new_n{:}),'Delimiter','tab');
																					    new_n=strcat(newStr(1),'-col_except_diag.txt') 
																					        writematrix(col,sprintf('%s',new_n{:}),'Delimiter','tab');
																						    new_n=strcat(newStr(1),'-values_except_diag.txt') 
																						        writematrix(v,sprintf('%s',new_n{:}),'Delimiter','tab');
																							    
																						end

