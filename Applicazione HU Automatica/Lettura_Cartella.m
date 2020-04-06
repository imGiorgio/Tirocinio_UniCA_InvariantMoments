%la lettura della cartella funziona, se non la trova restituisce un
%messaggio di errore
cartella = '/Users/macdigiorgio/Desktop/Dataset CT';
if ~isfolder(cartella)
  errorMessage = sprintf('Error: The following folder does not exist:\n%s', cartella);
  uiwait(warndlg(errorMessage));
  return;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('\n\n\n ------- Inizio ciclo -------\n\n');

%creo una lista di tutti i file all'interno della cartella con un
%determinato pattern
pattern = fullfile(cartella, '*.tiff');
immagini = dir(pattern);
for k = 1 : length(immagini)
   baseFileName = immagini(k).name;
   fullFileName = fullfile(cartella, baseFileName);
   fprintf(1, 'Ora sto leggendo %s\n', fullFileName);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Se messi all'interno del ciclo fanno vedere tutte e 39 le immagini in
%rapida sequenza direttamente a video

%imageArray = imread(fullFileName); 
%imshow(imageArray); %Mostro immagine a video
%drawnow;  %aggiorna lo schermo immediatamente
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
