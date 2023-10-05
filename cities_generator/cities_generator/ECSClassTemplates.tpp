template <typename TComponent>
TComponent* System::SingleCompFindQueue(std::string name)
{
    Queue* queue = GetQueueByName(name);
    return static_cast<TComponent*>((*queue->begin())[0]);
}
