import pickle
from protmapper import MappedSite

with open('output/pc_site_cache.pkl', 'rb') as f:
    cache = pickle.load(f)

new_cache = {}
for site_key, old_ms in cache.items():
    new_ms = MappedSite(up_id=old_ms.up_id,
                        error_code=None,
                        description=old_ms.description,
                        gene_name=old_ms.gene_name,
                        mapped_pos=old_ms.mapped_pos,
                        mapped_res=old_ms.mapped_res,
                        orig_pos=old_ms.orig_pos,
                        orig_res=old_ms.orig_res,
                        valid=old_ms.valid)
    new_cache[site_key] = new_ms

with open('output/pc_site_cache_new.pkl', 'wb') as f:
    pickle.dump(new_cache, f)
