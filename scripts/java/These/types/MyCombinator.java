package types;

import java.util.HashSet;
import java.util.Hashtable;
import java.util.Vector;

public class MyCombinator<U> {

	public Vector<HashSet<U>> liste = new Vector<HashSet<U>>();
	private Hashtable<U, Integer> dic = new Hashtable<U, Integer>();
	
	public void addLink(Vector<U> links) {
		if (links.size() == 0)
			return;
		int n = liste.size();
		HashSet<U> set = new HashSet<U>();
		for(U obj: links) {
			Integer j = dic.get(obj);
			if (j == null) {
				set.add(obj);
			} else if (j != n){
				HashSet<U> ens = liste.get(j);
				set.addAll(ens);
				for(U tmp: ens)
					dic.put(tmp, n);
				ens.clear();
			}
		}
		liste.add(set);
	}
	
}
