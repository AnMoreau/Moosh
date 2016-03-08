function eps=epsglass(lambda)
				%From "Fused Quartz Glass" range 346.69-3500nm
				% Rodney Spindler JOSA 44, 677 (1954)
if (lambda>346.69)&(lambda<3500)
eps=2.978645+0.008777808/(lambda^2*1e-6-0.010609)+84.06224/(lambda^2*1e-6-96);
else
disp('Erreur : hors du domaine de validité du modèle pour le verre')
end

end


