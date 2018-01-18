namespace MatchingAlgorithms{
  template <class TParticle, class TTrigger> // <3
  const std::tuple<TTrigger*, TParticle*, float, int>
  matchParticleWithL1Object
  (
    const TParticle &particle,
    const edm::Handle<std::vector<TParticle> > &particleCollectionHandle,
    const edm::Handle<BXVector<TTrigger> > &l1tObjectCollectionHandle,
    float dr2Max,
    bool performCrossCheck
  );

  template <class TParticle, class TTrigger> // <3
  const std::tuple<TTrigger *, TParticle *, float>
  matchL1ObjectWithParticle
  (
    const TTrigger &l1tObject,
    const edm::Handle<BXVector<TTrigger>> &l1tObjectCollectionHandle,
    const edm::Handle<std::vector<TParticle>> &particleCollectionHandle,
    float dr2Min,
    float dr2Max,
    bool performCrossCheck
  );
}