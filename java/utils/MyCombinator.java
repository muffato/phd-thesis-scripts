package utils;

import java.util.ArrayList;
import java.util.HashMap;

public class MyCombinator<U> {

	public ArrayList<ArrayList<U>> liste = new ArrayList<ArrayList<U>>();
	private HashMap<U, Integer> dic = new HashMap<U, Integer>();
	
	public void addElement(ArrayList<U> elt) {
		Integer n = new Integer(liste.size());
		for(U obj: elt) {
			dic.put(obj, n);
		}
		liste.add(elt);
	}
	
	public void addLink(U[] links) {
		if (links.length == 0)
			return;
		Integer n = new Integer(liste.size());
		ArrayList<U> set = new ArrayList<U>();
		for(U obj: links) {
			Integer j = dic.get(obj);
			if (j == null) {
				set.add(obj);
			} else if (j != n){
				ArrayList<U> ens = liste.get(j.intValue());
				set.addAll(ens);
				for(U tmp: ens)
					dic.put(tmp, n);
				ens.clear();
			}
		}
		liste.add(set);
	}
	
}
