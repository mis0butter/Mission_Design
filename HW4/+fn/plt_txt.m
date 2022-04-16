function plt_txt(txt)

pos = get(gca, 'position'); 

annotation('textbox', pos, ...
  'String', txt, ...
  'edgecolor', 'none');
axis off 

end 